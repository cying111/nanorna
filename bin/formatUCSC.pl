#!/usr/bin/perl -w

$input = shift @ARGV;

open $file, '<', "$input" or die "$0: Cannot open $file";
while ($line = <$file>){
  # chr 1-22 or X or Y
  if ($line =~ /^[0-9XY]/){
    $line =~ s/^/chr/;
  } elsif ($line =~ /^MT/){
    # mitochondrial
    $line =~ s/^MT/chrM/;
  }
  print "$line";
}
close $file;
