#!/usr/local/bin/perl

open(DATA, "SDF2asc $ARGV[0] vnorm vtan|");
while(<DATA>){
  ($vnorm, $vtan) = split;
  $nlines++;
  $sumnorm += $vnorm*$vnorm;
  $sumtan += $vtan*$vtan;
}

if ($nlines > 0){
  print "norm vel: ", sqrt($sumnorm/$nlines), 
	" tan vel: ", sqrt($sumtan/$nlines), "\n" ;
}else{
  print "no lines\n";
}
