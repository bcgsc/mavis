#!/projects/trans_scratch/software/perl/perl-5.20.3/bin/perl

package FusionProcessing;

#** @file
# Main file for processing Trans-ABySS results and connecting to ensembl
#*

use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';
use File::Basename;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
use Bio::EnsEMBL::ArchiveStableId;
use Try::Tiny;
use POSIX qw(strftime);
use TSV;
use JSON;


my $registry;
my $_version = '2.1.4';
my $_program = basename(__FILE__);
my $_install = dirname(abs_path(__FILE__));

main();

sub main
{
    my $outputfile;
    my $drug_target_file = '/projects/tumour_char/analysis_scripts/databases/processed_files/drug_target_tables/compiled_gene_drug_pathway.v1_2_4.tsv'; 
    my $best_transcript_file = 'ens69_best_transcript.txt';
    my $option_check = GetOptions(
        "output=s" => \$outputfile,
        "best_transcript_file" => \$best_transcript_file,
        "hugo_mapping_file" => \$drug_target_file
    );
    
    my $database_information =  {   -host => 'ensembl01.bcgsc.ca',
                                    -user => 'ensembl',
                                    -port => 3399,
                                    -pass => 'ensembl' };

    # set up the default filenames
    die "error: required argument --output not provided" if !defined $outputfile;

    $registry = 'Bio::EnsEMBL::Registry';
    $registry->load_registry_from_db(%$database_information);
    
    # read in the drug target file and generate a mapping for the ensembl gene id's
    my @required_column_names = ('ensid', 'hugo');
    print "loading: $drug_target_file\n";
    my ($header, $rows) = TSV::parse_input($drug_target_file, \@required_column_names);
    my %hugo_mapping = ();

    while (my $row = shift @$rows)
    {
        my $ensid = $row->{'ensid'};
        my $hugo = $row->{'hugo'};
        my @fields = split /;/, $hugo;
        $hugo_mapping{$ensid} = \@fields; 
    }
    
    @required_column_names = ('gene_id', 'transcript_id');
    print "loading: $best_transcript_file\n";
    ($header, $rows) = TSV::parse_input($best_transcript_file, \@required_column_names);
    my %best_transcript_mapping = ();

    while (my $row = shift @$rows)
    {
        my $ensid = $row->{'gene_id'};
        my $transcript = $row->{'transcript_id'};
        $best_transcript_mapping{$ensid} = $transcript;
    }
    
    
    
    # load all the different transcripts
    my $transcript_adaptor = $registry->get_adaptor('human', 'core', 'gene'); 
    my @glist = @{$transcript_adaptor->fetch_all()};
    my $counter = 1;
    my $total = scalar @glist;
    my $interval = $total / 100;
    print "loading $total genes\n";
    my %all_domains = ();
    my $time = localtime();
    my $jsons = {
        "hugo_mapping_file" => $drug_target_file,
        "best_transcript_file" => $best_transcript_file,
        "ensembl_version" => software_version(),
        "generation_time" => "$time",
        "script" => $_program,
        "script_version" => $_version,
        "genes" => []
    };

    while ( my $gene = shift @glist )
    {
        if ($counter % $interval == 0){
            my $percent = $counter * 100 / $total;
            print "$percent % complete\n";
        }
        
        my @tlist = @{$gene->get_all_Transcripts()};
        my $gid = $gene->stable_id();
        
        # get all hugo aliases for this ensembl gene
        my $hugo = "";
        if ( exists $hugo_mapping{$gid} ){
            $hugo = $hugo_mapping{$gid};
        }
        my $gjson = {
            "name" => $gid,
            "aliases" => $hugo,
            "transcripts" => [],
            "chr" => $gene->seq_region_name(),
            "start"=> $gene->start(),
            "end" => $gene->end(),
            "strand" => $gene->strand()
        };
        while ( my $t = shift @tlist )
        {
            my $tid = $t->stable_id();
            my $best = 'false';
            if ( exists $best_transcript_mapping{$gid} ){
                $best = 'true';
            }
            my $tjson = {
                "name" => $tid,
                "is_best_transcript" => $best,
                "exons" => [],
                "start" => $t->start(),
                "end" => $t->end(),
                "aliases" => [],
                "cdna_coding_start" => $t->cdna_coding_start(),
                "cdna_coding_end" => $t->cdna_coding_end(),
                "domains" => []
            };
            
            
            my $s_obj = $t->seq();
            my $s = $s_obj->seq();
            my $cds_start = $t->cdna_coding_start();
            my $cds_end = $t->cdna_coding_end();
            if ( !defined $cds_start || !defined $cds_end ){
                next;
            }
            # get all the refseq aliases for this ensembl transcript
            my @arr = @{$t->get_all_xrefs()};
            my @refseq = ();
            while (my $x  = shift @arr )
            {
                if ( ! ( $x->dbname() =~ /^RefSeq.*/ ) ){
                    next;
                }
                push(@{$tjson->{"aliases"}}, $x->display_id());
            }
            # get the translation start and end
            # get the domain coordinates (in amino-acids)
            # now add all of the domains
            my @domain_list = @{ $t->translation()->get_all_DomainFeatures() };
            my $domain_hash = {};
            for my $dom (@domain_list)
            {
                my $key = $dom->display_id(); # ensembl domain regions are split, group by display id
                my $curr = {
                    "regions" => [],
                    "desc" => $dom->idesc(),
                    "name" => $key
                };

                if(! exists $domain_hash->{$key})
                {
                    $domain_hash->{$key} = $curr;
                } else {
                    $curr = $domain_hash->{$key};
                }
                my $region = {"start" => $dom->start(), "end" => $dom->end() };
                push(@{$curr->{"regions"}},  $region);
            }
            foreach my $val (values %$domain_hash){
                push(@{$tjson->{"domains"}}, $val)
            }

            my @exon_list = @{ $t->get_all_Exons() };
            for my $ex (@exon_list)
            {
                my $exj = {
                    "start" => $ex->start(),
                    "end" => $ex->end(),
                    "name" => $ex->stable_id()
                };
                push(@{$tjson->{"exons"}}, $exj);
            }
            push(@{$gjson->{"transcripts"}}, $tjson);
        }
        if (scalar @{$gjson->{"transcripts"}} > 0){
            push(@{$jsons->{"genes"}}, $gjson);
            $counter = $counter + 1;
        }
    }
    open(my $fh, ">", $outputfile) or die "[ERROR] Could not open outputfile $outputfile\n";
    print "writeing: $outputfile\n";
    print $fh encode_json $jsons;
    close $fh;
    print "[$_program] [COMPLETE] status: Complete!\n";
}

