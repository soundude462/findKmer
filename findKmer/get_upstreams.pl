#!/usr/bin/perl

#This perl script extracts the most recent ensembl genes from a given database and creates a file of
#their upstream regions. chmod 777 get_upstreams.pl 

use strict;

use lib '/usr/local/bioperl/ensembl/modules';
use lib '/usr/local/bioperl/bioperl-live';
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Gene;



my($genefile, $contentsfile, $registry, @db_adaptors, $db_connection, $db_adaptor, $meta_container, $gene_adaptor, $len, $i, @stable_gene_ids, %db_match, @db_list, $upstream_size,
$gene, $other_name, $organism,@transcripts, $num_transcripts, $j, @exons, $num_exons, $exonfile, $db_version, $version, $chrom, $start, $end, $strand, $sqlfile, $out,
$sql_string, $release_num, $assembly, $genebuild, $transcriptsfile, $seq, $db, $userid, $connectionInfo, $dbh, $query, $dbh, $sth, %existing_releases, $seq_ref, %trans_genes,
    $ensembl_release, $org, $new_ensembl_id, $my_release_id, %contents_data, %gene_data, %transcript_data, $ensembl_key, @data, $transcript_key, $trans_ref, $seq, $outfile,
%existing_ensembl_versions, $latest_release, $ensembl_id, $transcript_key, $latest_ensembl, $seq_len, $cnt, @genes, $gene_ref, $ensembl_id, %gene_transcripts, $transcript,
%trans_info, $start_slice, $end_slice, $piece);

if (@ARGV >2 ) {
    $organism=$ARGV[0];
    $outfile=$ARGV[1];
    $upstream_size=$ARGV[2];

    $organism =~ s/\s/_/;




    $registry = 'Bio::EnsEMBL::Registry';

    $registry->load_registry_from_db(
        -host => 'ensembldb.ensembl.org',
        -user => 'anonymous'
    );

    $db_adaptor = $registry->get_DBAdaptor($organism, "core");

    my @db_adaptors = @{ $registry->get_all_DBAdaptors() };

    $new_ensembl_id="";

    foreach my $db_adaptor (@db_adaptors) {
        my $db_connection = $db_adaptor->dbc();

        #printf(
        #"species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
        #$db_adaptor->species(),   $db_adaptor->group(),
        #$db_connection->dbname(), $db_connection->host(),
        #$db_connection->port()
        #);

        
        if (($db_adaptor->species() =~ /$organism/) && ($db_adaptor->group() eq "core")) {
            $new_ensembl_id = $db_connection->dbname();
            $new_ensembl_id =~ s/^$organism//;
            $new_ensembl_id =~ s/^_core_//;
            $new_ensembl_id =~ s/_\d{1,3}$//;
            
            print "Found New Ensembl ID $new_ensembl_id for organism $organism\n";
        }
        
      
    }

    if ($new_ensembl_id eq "") {
        print "ERROR: could not find an ensembl release for organism $organism\n";
    }
    else {
        open(WRITEDATA, ">" . $outfile) or die;
        my $trans_adaptor  = $registry->get_adaptor( $organism, 'Core', 'Transcript' );
        @stable_gene_ids = @{$trans_adaptor->list_stable_ids()};
        $trans_ref = $trans_adaptor->fetch_all_by_stable_id_list	(\@stable_gene_ids);
    
      
        while ($transcript = shift @{$trans_ref} ) {
            $ensembl_id=$transcript->display_id();
            
            
            if ($transcript->biotype() eq "protein_coding") {
                $version=$transcript->version();
                $seq=$transcript->slice();
                if ($seq) {
                    $end=$transcript->start();
                    $start=$end-$upstream_size;
                    $piece = $seq->subseq($start, $end);
                    
                    print WRITEDATA ">", $ensembl_id, "\n", $piece, "\n";
                }
            }
  
        
        }
        close(WRITEDATA);
    }
}
else {
    print "Usage: get_upstreams.pl organism_name <output file> upstream_length\n";
}
