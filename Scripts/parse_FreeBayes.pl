use strict;
use warnings;



open (READ, "<$ARGV[0]"); #VCF;

while (<READ>) {
	my $line = $_;
	if ($line =~ /^#/) {print $line;next}
	chomp;
	my @mas = split/\t/;
	my @format = split/:/, $mas[8];
	my @field  = split/:/, $mas[9];
	my $GT;
	for (my $i = 0; $i < scalar @format; $i++) {
		if ($format[$i] eq 'GT') {
			$GT = $i;
			}
		}
	next if $field[$GT] eq '0/0';
	next if $mas[5] < 5;
	$mas[6] = 'PASS';
	$line = join("\t", @mas);
	print "$line\n";
	next;
	}

close READ;
