import re

from pysam import VariantFile

from ..constants import COLUMNS, ORIENT, SVTYPE
from ..util import DEVNULL
from ..interval import Interval
from .constants import SUPPORTED_TOOL


def parse_bnd_alt(alt):
    """
    parses the alt statement from vcf files using the specification in vcf 4.2/4.2.

    Assumes that the reference base is always the outermost base (this is based on the spec and also manta results as
    the spec was missing some cases)

    r = reference base/seq
    u = untemplated sequence/alternate sequence
    p = chromosome:position

    | alt format   | orients |
    | ------------ | ------- |
    | ru[p[        | LR      |
    | [p[ur        | RR      |
    | ]p]ur        | RL      |
    | ru]p]        | LL      |
    | u.           |         |
    | .u           |         |
    """
    # ru[p[
    match = re.match(r"^(?P<ref>\w)(?P<useq>\w*)\[(?P<chr>[^:]+):(?P<pos>\d+)\[$", alt)
    if match:
        return (
            match.group("chr"),
            int(match.group("pos")),
            ORIENT.LEFT,
            ORIENT.RIGHT,
            match.group("ref"),
            match.group("useq"),
        )
    # [p[ur
    match = re.match(r"^\[(?P<chr>[^:]+):(?P<pos>\d+)\[(?P<useq>\w*)(?P<ref>\w)$", alt)
    if match:
        return (
            match.group("chr"),
            int(match.group("pos")),
            ORIENT.RIGHT,
            ORIENT.RIGHT,
            match.group("ref"),
            match.group("useq"),
        )
    # ]p]ur
    match = re.match(r"^\](?P<chr>[^:]+):(?P<pos>\d+)\](?P<useq>\w*)(?P<ref>\w)$", alt)
    if match:
        return (
            match.group("chr"),
            int(match.group("pos")),
            ORIENT.RIGHT,
            ORIENT.LEFT,
            match.group("ref"),
            match.group("useq"),
        )
    # ru]p]
    match = re.match(r"^(?P<ref>\w)(?P<useq>\w*)\](?P<chr>[^:]+):(?P<pos>\d+)\]$", alt)
    if match:
        return (
            match.group("chr"),
            int(match.group("pos")),
            ORIENT.LEFT,
            ORIENT.LEFT,
            match.group("ref"),
            match.group("useq"),
        )
    else:
        raise NotImplementedError("alt specification in unexpected format", alt)


def convert_record(record, record_mapping={}, log=DEVNULL):
    """
    converts a vcf record

    Note:
        CT = connection type, If given this field will be used in determining the orientation at the breakpoints.
        From https://groups.google.com/forum/#!topic/delly-users/6Mq2juBraRY, we can expect certain CT types for
        certain event types
            - translocation/inverted translocation: 3to3, 3to5, 5to3, 5to5
            - inversion: 3to3, 5to5
            - deletion: 3to5
            - duplication: 5to3
    """
    records = []
    for alt in record.alts if record.alts else [None]:
        info = {}
        for key in record.info.keys():
            try:
                value = record.info[key]
            except UnicodeDecodeError as err:
                log("Ignoring invalid INFO field {} with error: {}".format(key, err))
            else:
                try:
                    value = value[0] if len(value) == 1 else value
                except TypeError:
                    pass  # anything non-tuple
            info[key] = value

        std_row = {}
        if record.id and record.id != "N":  # to account for NovoBreak N in the ID field
            std_row["id"] = record.id

        if info.get("SVTYPE", None) == "BND":
            if "." in alt:
                # ru.
                match = re.match(r"^(?P<ref>\w)(?P<useq>\w*)\.$", alt)
                if match:
                    chr2 = record.chrom
                    std_row[COLUMNS.break1_orientation] = ORIENT.NS
                    std_row[COLUMNS.break2_orientation] = ORIENT.NS
                    std_row[COLUMNS.untemplated_seq] = match.group("useq")
                # .ur
                match = re.match(r"^\.(?P<useq>\w*)(?P<ref>\w)$", alt)
                if match:
                    chr2 = record.chrom
                    std_row[COLUMNS.break1_orientation] = ORIENT.NS
                    std_row[COLUMNS.break2_orientation] = ORIENT.NS
                    std_row[COLUMNS.untemplated_seq] = match.group("useq")
                end = record.stop

            else:
                chr2, end, orient1, orient2, ref, alt = parse_bnd_alt(alt)
                std_row[COLUMNS.break1_orientation] = orient1
                std_row[COLUMNS.break2_orientation] = orient2
                std_row[COLUMNS.untemplated_seq] = alt
                if record.ref != ref:
                    raise AssertionError(
                        "Expected the ref specification in the vcf record to match the sequence "
                        "in the alt string: {} vs {}".format(record.ref, ref)
                    )
        else:
            chr2 = info.get("CHR2", record.chrom)
            end = record.stop
            if (
                alt
                and record.ref
                and re.match(r"^[A-Z]+$", alt)
                and re.match(r"^[A-Z]+", record.ref)
            ):
                std_row[COLUMNS.untemplated_seq] = alt[1:]
                size = len(alt) - len(record.ref)
                if size > 0:
                    std_row[COLUMNS.event_type] = SVTYPE.INS
                elif size < 0:
                    std_row[COLUMNS.event_type] = SVTYPE.DEL
        std_row.update({COLUMNS.break1_chromosome: record.chrom, COLUMNS.break2_chromosome: chr2})
        if info.get(
            "PRECISE", False
        ):  # DELLY CI only apply when split reads were not used to refine the breakpoint which is then flagged
            std_row.update(
                {
                    COLUMNS.break1_position_start: record.pos,
                    COLUMNS.break1_position_end: record.pos,
                    COLUMNS.break2_position_start: end,
                    COLUMNS.break2_position_end: end,
                }
            )
        else:
            convert_imprecise_breakend(std_row, record, info, end)

        if "SVTYPE" in info:
            std_row[COLUMNS.event_type] = info["SVTYPE"]

        try:
            orient1, orient2 = info["CT"].split("to")
            connection_type = {"3": ORIENT.LEFT, "5": ORIENT.RIGHT, "N": ORIENT.NS}
            std_row[COLUMNS.break1_orientation] = connection_type[orient1]
            std_row[COLUMNS.break2_orientation] = connection_type[orient2]
        except KeyError:
            pass
        std_row.update(
            {k: v for k, v in info.items() if k not in {"CHR2", "SVTYPE", "CIPOS", "CIEND", "CT"}}
        )
        records.append(std_row)
    return records


def convert_imprecise_breakend(std_row, record, info, bp_end):
    """
    Handles IMPRECISE calls, that leveraged uncertainty from the CIPOS/CIEND/CILEN fields.

    bp1_s = breakpoint1 start
    bp1_e = breakpoint1 end
    bp2_s = breakpoint2 start
    bp2_e = breakpoint2 end

    Insertion and deletion edge case - in which bp1_e > bp2_s
    E.g bp1_s = 1890, bp1_e = 2000, bp2_s = 1900, bp2_e = 1900.
    break1 ------------------------=======================--------------
    break2 ------------------------==========---------------------------

    Insertion edge case - in which bp1_e > bp1_s
    E.g bp1_s = 1890, bp1_e = 1800, bp2_s = 1800, bp2_e = 1800.
    break1 ------------------------==-----------------------------------
    break2 ------------------------=------------------------------------

    Insertion edge case - in which bp1_s > bp2_s
    E.g bp1_s = 1950, bp1_e = 2000, bp2_s = 1900, bp2_e = 3000.
    break1 ------------------------==-----------------------------------
    break2 -----------------------========------------------------------
    """
    is_intrachromosomal = std_row["break1_chromosome"] == std_row["break2_chromosome"]

    std_row.update(
        {
            COLUMNS.break1_position_start: max(1, record.pos + info.get("CIPOS", (0, 0))[0]),
            COLUMNS.break1_position_end: record.pos + info.get("CIPOS", (0, 0))[1],
            COLUMNS.break2_position_start: max(1, bp_end + info.get("CIEND", (0, 0))[0]),
            COLUMNS.break2_position_end: bp_end + info.get("CIEND", (0, 0))[1],
        }
    )

    if is_intrachromosomal and (
        Interval.overlaps(
            (std_row["break1_position_start"], std_row["break1_position_end"]),
            (std_row["break2_position_start"], std_row["break2_position_end"]),
        )
    ):
        if std_row.get("event_type") != "BND":
            std_row["break2_position_start"] = max(
                std_row["break1_position_start"], std_row["break2_position_start"]
            )
            std_row["break1_position_end"] = min(
                std_row["break1_position_end"], std_row["break2_position_end"]
            )
            std_row["break1_position_start"] = min(
                std_row["break1_position_start"], std_row["break2_position_end"]
            )
    if std_row["break1_position_end"] == 0 and std_row["break1_position_start"] == 1:
        # addresses cases where pos = 0 and telomeric BND alt syntax https://github.com/bcgsc/mavis/issues/294
        std_row.update({"break1_position_end": 1})
    if std_row["break2_position_end"] == 0 and std_row["break2_position_start"] == 1:
        std_row.update({"break2_position_end": 1})

    if is_intrachromosomal and (
        std_row["break2_position_start"] > std_row["break2_position_end"]
        or std_row["break1_position_start"] > std_row["break1_position_end"]
    ):
        if "event_type" in std_row and std_row["event_type"] != "BND":
            raise ValueError(
                f'Improper entry. One of the following breakpoints start is greater than breakpoint end: Breakpoint1_start: {std_row["break1_position_start"]}, Breakpoint1_end: {std_row["break1_position_end"]} Breakpoint2_start: {std_row["break2_position_start"]}, Breakpoint2_end: {std_row["break2_position_end"]} This call has been dropped.'
            )

    if (
        None not in (record.pos, record.info.get("END"))
        and record.pos > record.info.get("END")
        and is_intrachromosomal
    ):
        raise ValueError(
            f'Improper entry. Starting position ({record.pos}) cannot be greater than ending position ({record.info.get("END")}). This call has been dropped.'
        )


def convert_file(input_file: str, file_type: str, log):
    """process a VCF file

    Args:
        input_file: the input file name
        file_type: the input type

    Raises:
        err: [description]
    """
    rows = []
    vfile = VariantFile(input_file)
    try:
        vfile.header.info.add("END", number=1, type="Integer", description="End of the interval")
    except ValueError:
        pass

    for vcf_record in vfile.fetch():
        try:
            rows.extend(convert_record(vcf_record, log=log))
        except Exception as err:
            if file_type not in [SUPPORTED_TOOL.STRELKA, SUPPORTED_TOOL.MUTECT]:
                raise err
            else:
                log("Ignoring", vcf_record)
    return rows
