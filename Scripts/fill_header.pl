use strict;
use warnings;
use Dir::Self;

my $header_path = __DIR__ . "/header.vcf";
open (INFO, '<', "$header_path") or die;
my $info = '';
while (<INFO>) {
	$info = $info.$_;
}
close INFO;

my $k = 0;
while (<STDIN>) {
	if ((/INFO=/) && ($k == 0)) {
		print $info;
		print;
		$k = 1;
	} else {
		print;
	}
}

