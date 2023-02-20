$panel = $ARGV[0];
chomp $panel;
open (BED, '<', $panel) or die;
%positions = ();
%amplicons = ();
%lengths = ();

while (<BED>) {
	next if (/^track/);
	@a = split/\t/;
	next if ($a[0] eq '');
	$chr = $a[0];
	$left = $a[1];
	$right = $a[2];
	$amplicons{$chr."\t".$left."\t".$right} = 0;
	$lengths{$chr."\t".$left."\t".$right} = $right-$left+1;
	for ($i=$left;$i<=$right;$i++) {
		$positions{$chr."\t".$i}{"$chr"."\t"."$left"."\t"."$right"} = 0;
		#print "$chr"."\t"."$i"."\n";
	}
}
close BED;
#print "ANALYSIS\n";
while (<STDIN>) {
	chomp;
	@b = split/\t/;
	$ampl = '';
	if (exists($positions{$b[0]."\t".$b[1]})) {
		foreach $p (keys %{$positions{$b[0]."\t".$b[1]}}) {
			#print;
			#print "\n";
			$amplicons{$p} += $b[2];
			#$lengths{$ampl} += 1;
		}
	}
}

foreach $k (sort keys %amplicons) {
	if ($lengths{$k} == 0) {
		print "$k\tYaba-daba-doo\t0\n";
	} else {
		$avg = $amplicons{$k}/$lengths{$k};
		print "$k\tYaba-daba-doo\t$avg\n";
	}
}

