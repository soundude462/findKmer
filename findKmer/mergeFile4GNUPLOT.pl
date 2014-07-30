
#!/usr/bin/perl

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
