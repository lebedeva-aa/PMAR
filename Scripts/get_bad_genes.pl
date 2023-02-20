use Dir::Self;
$current_dir = __DIR__;
$size = $ARGV[0];
$uniformity = $ARGV[1];
$table = $ARGV[2];

open (BED, '<', "$current_dir/../conf/WES.bed") or die;
%bed = ();
while (<BED>) {
	chomp;
	@b = split/\t/;
	$name = $b[0]."\t".$b[1]."\t".$b[2];
	#@c = split/;/, $b[5];
	#@d = split/=/, $c[0];
	#$bed{$name} = $d[1];
	$bed{$name} = $b[3];
}
close BED;

open (SIZE, '<', $size) or die;
while (<SIZE>) {
	print "# Total reads:\t$_";
}
close SIZE;

open (UNI, '<', $uniformity) or die;
@file = <UNI>;
print "# Bases in targeted regions:\t$file[0]";
print "# Uniformity:\t$file[1]\n";
close UNI;

open (IN, '<', $table) or die;
while (<IN>) {
	print if (/^#/);
	next if (/^#/);
	chomp;
	$l = $_;
	@a = split/\t/, $l;
	$n = $a[0]."\t".$a[1]."\t".$a[2];
	$gene = $bed{$n};
	print "$l\t$gene\n";
}
close IN;
