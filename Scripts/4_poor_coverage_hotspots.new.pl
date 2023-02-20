$depth = $ARGV[0];

open (LIST, '<', '/home/aod/AOD/Clinical_exome/Scripts/exome.hotspots.txt') or die;
%list = ();
while (<LIST>) {
	chomp;
	@q = split/\t/;
	$list{"$q[0]\t$q[1]"} = "$q[2]\t$q[3]";
}
close LIST;

open (DEP, '<', $depth) or die;
%indels = ();
print "Chromosome\tPosition\tGene\tAmino acid substitution\tDepth of coverage\n";
while (<DEP>) {
	chomp;
	@a = split/\t/;
	$r = $list{"$a[0]\t$a[1]"};
	print "$a[0]\t$a[1]\t$r\t$a[2]\n" if ((exists($list{"$a[0]\t$a[1]"})) && ($a[2] <= 30));
}
close DEP;

