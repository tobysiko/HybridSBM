#!/usr/bin/perl

@lines=`grep -v '#' n0/temperature.info`;
for ($i=0;$i<@lines;++$i) {
    $line=$lines[$i];
    @parts=split(/\s+/,$line);
    $beta[$i]=$parts[3];
    $temp[$i]=$parts[1];
    $tempk[$i]=$parts[2];
#    print "$i  $temp[$i]\n";
}

$file_template="averages";
if (@ARGV ==1) {
    $file_template=$ARGV[0];
} 

@files=<n*/$file_template>;
 
open(fp,"<n0/$file_template");
@lines=<fp>;
$obsno=0;
$linecount=0;
foreach $line (@lines) {
    
    ($f1,$f2,$f3) = split(/\s+/,$line);
    if ($f1=~m/\#/ && $linecount>=2) {
	@obsname[$obsno]=$f2;
	@obsstart[$obsno]=$linecount+1;
	if ($obsno!=0) {
	    @obsend[$obsno-1]=$linecount;
	}
	$obsno++;
    }
    $linecount++;
}
$obsend[$obsno-1]=$linecount;
close(fp);
$obscount=0;
foreach $obs (@obsname) {
    $stline=$obsstart[$obscount];
    $ndline=$obsend[$obscount];
    $ifile=0;
    foreach $file (@files) {
	open(fp,"<$file");
	@lines=<fp>;
	for ($i=$stline;$i<$ndline;++$i) {
	    $line = @lines[$i];
	    ($tindex,$mean,$stddev)=split(/\s/,$line);
#	print "$tindex $mean\n";
	    if ($ifile==0) {
		@mean_mean[$tindex]=$mean;
		@mean_stddev[$tindex]=$stddev;
		@stddev_mean[$tindex]=$mean*$mean;
	    } else {
		@mean_mean[$tindex]=$mean_mean[$tindex]+$mean;
		@mean_stddev[$tindex]=$mean_stddev[$tindex]+$stddev;
		@stddev_mean[$tindex]=$stddev_mean[$tindex]+$mean*$mean;
	    }
	}
	$ifile++;
	close(fp);
    }
    open(av,">$obs\.avg");
    print "Observable: $obs starts: $stline ends: $ndline\n";
    print av "# Index T(Kelvin) T(model)  Beta  mean(mean)  mean(std.dev.) std.dev.(mean) \n";
    for ($i=0;$i<=$#mean_mean;++$i) {
	$mn=$mean_mean[$i]=$mean_mean[$i]/$ifile;
	$mean_stddev[$i]=$mean_stddev[$i]/$ifile;
	$stddev_mean[$i]=sqrt($stddev_mean[$i]/$ifile-$mn*$mn);
	print av "$i  $tempk[$i]   $temp[$i]   $beta[$i]   $mean_mean[$i]   $mean_stddev[$i]   $stddev_mean[$i]\n";
	print "$i  $tempk[$i]  $temp[$i]   $beta[$i]   $mean_mean[$i]   $mean_stddev[$i]   $stddev_mean[$i]\n";
    }
    close(av);
    $obscount++;
}


