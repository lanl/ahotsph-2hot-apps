#!/usr/local/bin/perl

open(DATA, "SDFcvt -bf $ARGV[0] vnorm vtan|");
$maxnorm = 0.0;
$maxtan = 0.0;
while(read(DATA, $buf, 8)){
  ($vnorm, $vtan) = unpack("ff", $buf);
  $nlines++;
  $sumnorm += $vnorm*$vnorm;
  $sumtan += $vtan*$vtan;
  if( $vnorm > $maxnorm ){
      $maxnorm = $vnorm;
  }
  if( $vtan > $maxtan ){
      $maxtan = $vtan;
  }
}

if ($nlines > 0){
  print "norm vel rms: ", sqrt($sumnorm/$nlines), " max: ", $maxnorm, "\n", 
	" tan vel rms: ", sqrt($sumtan/$nlines), " max: ", $maxtan, "\n" ;
}else{
  print "no lines\n";
}
