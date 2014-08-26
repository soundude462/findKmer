
#!/usr/bin/perl

#This perl script reads two files, one file has more lines than the other.
#Both files have a DNA sequence as the first thing on each line
#The sequence is seperated by a comma
#This perl script checks to see if the sequence on each line is the same line for line in each file
#If the sequence does not match then it drops that sequence from the genome file and drops to the next line of the genome file
#This continues until a match is found between the genome file and the upstream file.
#We do this because there are sequences that occur in the genome that do not occur in the upstream.
#Our research is only interested in what occurs in the upstream region so genome only sequences are irrelevant. 
 
#chmod 777 get_upstreams.pl 

use strict;
my($file1, $file2, %kmers1, %allkmers, $line, @data, $outfile, $key);

$file1 =$ARGV[0];
$file2 =$ARGV[1];
$outfile =$ARGV[3];

open(READDATA, "<". $file1) or die;

while($line = <READDATA>) {
$line =~ s/\n//;
@data=split(/\t/. $line);
$kmers1{$data[0]}= $data[0] . "\t". $data[1] . "\t" . $data[2];
}
close(READDATA);

open(READDATA, "<". $file2) or die;

while($line = <READDATA>) {
$line =~ s/\n//;
@data=split(/\t/. $line);
if (exists $kmers1{$data[0]}) {
$allkmers{$data[0]}= $kmers1{$data[0]} .  "\t" . $data[2];
}
}
close(READDATA);

open (WRITEDATA, ">" . $outfile) or die;
foreach $key (keys %allkmers) {print WRITEDATA $allkmers{$key}, "\n";}
close(WRITEDATA);
