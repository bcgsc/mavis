import re

from ..constants import STRAND
from ..interval import Interval


class ReferenceName(str):
    """
    Class for reference sequence names. Ensures that hg19/hg38 chromosome names match.

    Example:
        >>> ReferenceName('chr1') == ReferenceName('1')
        True
    """
    def __eq__(self, other):
        options = {str(self)}
        if self.startswith('chr'):
            options.add(str(self[3:]))
        else:
            options.add('chr' + str(self))
        return other in options

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(re.sub('^chr', '', str(self)))

    def __lt__(self, other):
        self_std_repr = self if not self.startswith('chr') else self[3:]
        other_std_repr = other if not other.startswith('chr') else other[3:]
        return str.__lt__(self_std_repr, other_std_repr)

    def __gt__(self, other):
        self_std_repr = self if not self.startswith('chr') else self[3:]
        other_std_repr = other if not other.startswith('chr') else other[3:]
        return str.__gt__(self_std_repr, other_std_repr)

    def __ge__(self, other):
        if self == other:
            return True
        return self.__gt__(other)

    def __le__(self, other):
        if self == other:
            return True
        return self.__lt__(other)


class BioInterval:

    def __init__(self, reference_object, start, end=None, name=None, seq=None, data=None, strand=None):
        """
        Args:
            reference_object: the object this interval is on
            start (int) start of the interval (inclusive)
            end (int): end of the interval (inclusive)
            name: optional
            seq (str): the seq relating to this interval

        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> b[0]
            12572784
            >>> b[1]
            12578898
        """
        start = int(start)
        end = int(end) if end is not None else None
        data = {} if data is None else data
        self.reference_object = reference_object
        self.name = name
        self.position = Interval(start, end, number_type=int)
        self.seq = seq if not seq else str(seq.upper())
        self.data = {}
        self.data.update(data)
        self.strand = strand

    @property
    def start(self):
        """*int*: the start position"""
        return self.position.start

    @property
    def end(self):
        """*int*: the end position"""
        return self.position.end

    def __getitem__(self, index):
        return Interval.__getitem__(self, index)

    def __len__(self):
        """
        Example:
            >>> b = BioInterval('1', 12572784, 12578898, 'q22.2')
            >>> len(b)
            6115
        """
        return self.position.length()

    def key(self):
        """:class:`tuple`: a tuple representing the items expected to be unique. for hashing and comparing"""
        return (self.reference_object, self.position, self.seq, self.name)

    def __eq__(self, other):
        if not hasattr(other, 'key'):
            return False
        return self.key() == other.key()

    def __lt__(self, other):
        if other.reference_object and self.reference_object and other.reference_object != self.reference_object:
            if self.reference_object < other.reference_object:
                return True
            return False
        elif self.position < other.position:
            return True
        return False

    def __hash__(self):
        return hash(self.key())

    def get_seq(self, reference_genome=None, ignore_cache=False):
        """
        get the sequence for the current annotation object

        Raises:
            NotImplementedError: abstract method
        """
        raise NotImplementedError('abstract method must be overidden')

    def get_strand(self):
        """
        pulls strand information from the current object, or follows reference
        objects until the strand is found

        Returns:
            STRAND: the strand of this or any of its reference objects

        Raises:
            AttributeError: raised if the strand is not set on this or any of its reference objects
        """
        try:
            if self.strand is not None:
                return self.strand
        except AttributeError:
            pass
        tried = set()
        parent = self.reference_object
        while True:
            if parent in tried:
                break
            try:
                if parent.strand is not None:
                    return parent.strand
            except AttributeError:
                pass
            tried.add(parent)
            try:
                parent = parent.reference_object
            except AttributeError:
                break
        raise AttributeError('strand has not been defined', self, self.strand, tried)

    @property
    def is_reverse(self):
        """True if the gene is on the reverse/negative strand.

        Raises:
            AttributeError: if the strand is not specified
        """
        if self.get_strand() == STRAND.NEG:
            return True
        elif self.get_strand() == STRAND.POS:
            return False
        else:
            raise AttributeError('strand has not been defined', self)

    def get_chr(self):
        """
        pulls chromosome information from the current object, or follows reference
        objects until the chromosome is found

        Returns:
            str: the chromosome of this or any of its reference objects

        Raises:
            AttributeError: raised if the chromosome is not set on this or any of its reference objects
        """
        try:
            if self.chr is not None:
                return self.chr
        except AttributeError:
            pass

        tried = set()
        parent = self.reference_object
        while True:
            if parent in tried:
                break
            try:
                if parent.chr is not None:
                    return parent.chr
            except AttributeError:
                pass
            tried.add(parent)
            try:
                parent = parent.reference_object
            except AttributeError:
                break
        raise AttributeError('chr has not been defined')

    def to_dict(self):
        """
        creates a dictionary representing the current object

        Returns:
            :class:`dict` by :class:`str`: the dictionary of attribute values
        """
        dict_result = {
            'name': self.name,
            'start': self.start,
            'end': self.end,
            'type': self.__class__.__name__,
            'seq': self.seq,
            'data': self.data
        }
        try:
            dict_result['reference_object'] = self.reference_object.name
        except AttributeError:
            dict_result['reference_object'] = str(self.reference_object)
        return dict_result

    def __repr__(self):
        cls = self.__class__.__name__
        refname = self.reference_object
        try:
            refname = self.reference_object.name
        except AttributeError:
            pass
        return '{}({}:{}-{}, name={})'.format(cls, refname, self.start, self.end, self.name)
