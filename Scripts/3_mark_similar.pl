@a = <STDIN>;
%hash = ();
foreach $t (@a) {
	next if ($t =~ /^Consequence/);
	@q = split/\t/, $t;
	($chr,$p) = split/:/, $q[3];
	$p =~ /\d+/;
	$pos = $&;
	#print "$pos\n";
	$hash{"$chr\t$pos"} += 1;
}

foreach $t (@a) {
	$before = '';
	$after = '';
	@q = split/\t/, $t;
	if ($t =~ /^Consequence/) {
		for ($y=0;$y<=23;$y++){
			$before = $before.$q[$y]."\t";
		}
		print "$before";
		print "Number of mutations in position";
		for ($y=24;$y<=$#q;$y++){
			$after = $after."\t".$q[$y];
		}
		print "$after";
		next;
	}
	for ($y=0;$y<=23;$y++){
		$before = $before.$q[$y]."\t";
	}
	print "$before";
	($chr,$p) = split/:/, $q[3];
	$p =~ /\d+/;
	$pos = $&;
	$l = $hash{"$chr\t$pos"};
	print "$l";
	for ($y=24;$y<=$#q;$y++){
		$after = $after."\t".$q[$y];
	}
	print "$after";
}
