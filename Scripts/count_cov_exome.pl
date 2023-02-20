$total = 0;
$zero = 0;
$one = 0;
$ten = 0;
$twenty = 0;
$hundred = 0;
$fivehu = 0;
@bad = ();
$depth = 0;
while (<STDIN>) {
	@a = split/\t/;
	next if ($a[2]-$a[1] <= 10);
	$zero += 1 if ($a[4] == 0);
	$one += 1 if ($a[4] < 1);
	$ten += 1 if ($a[4] < 10);
	$twenty += 1 if ($a[4] < 20);
	$p_30 += 1 if ($a[4] < 30);
	$p_50 += 1 if ($a[4] < 50);
	$hundred += 1 if ($a[4] < 100);
	$fivehu += 1 if ($a[4] < 500);
	push (@bad, $_) if ($a[4] < 20);
	$depth += $a[4];
	$total += 1;
}
$avg = $depth/$total;
$perc_ze = (1-$zero/$total)*100;
$perc_on = (1-$one/$total)*100;
$perc_te = (1-$ten/$total)*100;
$perc_tw = (1-$twenty/$total)*100;
$perc_hu = (1-$hundred/$total)*100;
$perc_fh = (1-$fivehu/$total)*100;
$p_30 = (1-$p_30/$total)*100;
$p_50 = (1-$p_50/$total)*100;
print "# Total regions: $total\n# Average depth: $avg\n# Not covered: $zero\n# Less x1: $one\n# Less x10: $ten\n# Less x20: $twenty\n# Less x100: $hundred\n# Less x500: $fivehu\n# Percentage of covered: $perc_ze\n# Percentage of higher x1: $perc_on\n# Percentage of higher x10: $perc_te\n# Percentage of higher x20: $perc_tw\n#Percentage of higher x30: $p_30\n#Percentage of higher x50: $p_50\n# Percentage of higher x100: $perc_hu\n# Percentage of higher x500: $perc_fh\n";
foreach $n (@bad) {
	print $n;
}
