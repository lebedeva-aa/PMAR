$panel = $ARGV[0];
chomp $panel;
open (BED, '<', $panel) or die;
%positions = ();

while (<BED>) {
	next if (/^track/);
	@a = split/\t/;
	next if ($a[0] eq '');
	chomp;
	$chr = $a[0];
	$left = $a[1];
	$right = $a[2];
	$percentile = $a[4]*0.2;
	for ($i=$left;$i<=$right;$i++) {
		$positions{$chr."\t".$i} = $percentile;
		#print "$chr"."\t"."$i"."\n";
	}
}
close BED;
@pos = keys %positions;
$all = @pos;
$sum = 0;
while (<STDIN>) {
	chomp;
	@b = split/\t/;
	$sum += 1 if (exists($positions{"$b[0]\t$b[1]"}) && ($b[2] >= $positions{"$b[0]\t$b[1]"}));
}
$uniformity = $sum/$all*100;
print "$all\n$uniformity\n";
