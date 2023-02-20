#%hash = ();
while (<STDIN>) {
	if (/^#/) {
		print;
	} else {
		@a = split/\t/;
		if ($a[6] eq 'PASS'){
			chomp $a[9];
			@c = split/:/, $a[8];
			@b = split/:/, $a[9];
			$k = 0;
			$n = 0;
			$chr = $a[0].':'.$a[1];
			$chr =~ s/chr//;
			$af = VCFinfo($a[7], 'AF') || 0;
			$dp = VCFinfo($a[7], 'DP') || 0;
			$ao = VCFinfo($a[7], 'AO') || 0;
			$saf = VCFinfo($a[7], 'SAF') || 0;
			$sar = VCFinfo($a[7], 'SAR') || 0;
			$srf = VCFinfo($a[7], 'SRF') || 0;
			$srr = VCFinfo($a[7], 'SRR') || 0;
			for ($i=0;$i<=$#c;$i++){
				if ($c[$i] eq 'AD'){
					@d = split/,/, $b[$i];
					$ao = $d[1];
					$af = $d[1]/($d[1]+$d[0]);
				}
				if ($c[$i] eq 'DP'){
					$dp = $b[$i];
				} elsif ($c[$i] eq 'DPI'){
					$dp = $b[$i];
				}
				if ($c[$i] eq 'ADF'){
					@d = split/,/, $b[$i];
					$srf = $d[0];
					$saf = $d[1];
				}
				if ($c[$i] eq 'ADR'){
					@d = split/,/, $b[$i];
					$srr = $d[0];
					$sar = $d[1];
				}
			}
			print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\tAF=$af;DP=$dp;AO=$ao;SAF=$saf;SAR=$sar;SRF=$srf;SRR=$srr\t$a[8]\t$a[9]\n";
		}
	}
}


sub VCFinfo {
        my $line = shift;
        my $field = shift;
        my @info = split/;/, $line;
        foreach my $arg (@info) {
                if ($arg =~ /^$field=(\S+)$/) {
                        return $1;
                        }
                }
        return undef;
        }

