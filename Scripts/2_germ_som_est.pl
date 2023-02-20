
#Используется только после добавления в начало 'AF_	DP	DB	'!!! Все названия столбцов должны соответствовать содержанию!!!

$af_i = 0;
$dp_i = 0;
$ex_v_i = 0;


while (<STDIN>){
	if ($_ !~ /^Consequence/) {
		@a = split/\t/;
		$before = '';
		$after = '';
		for ($y=0;$y<=23;$y++){
			$before = $before.$a[$y]."\t";
		}
		print "$before";
		$af = $a[$af_i];
		$cosm = $a[$cosm_i];
		$rs =  $a[$rs_i];
		$cov = $a[$dp_i]*$a[$af_i];
		if (($af >= 0.4) && ($af <= 0.6) && ($rs =~ /rs/) && ($cosm !~ /COSM/) && ($cov >= 5)) {
			print "Germline";
		} elsif (($af >= 0.85) && ($rs =~ /rs/) && ($cosm !~ /COSM/) && ($cov >= 5)) {
			print "Germline";
		} elsif (($af >= 0.4) && ($af <= 0.60) && ($cosm =~ /COSM/) && ($cov >= 5)) {
			print "Possibly somatic";
		} elsif (($af >= 0.85)  && ($cosm =~ /COSM/) && ($cov >= 5)) {
			print "Possibly somatic";
		} elsif (($af < 0.4)  && ($cosm =~ /COSM/) && ($cov >= 5)) {
			print "Somatic";
		} elsif (($af < 0.4)  && ($rs =~ /rs/) && ($cosm !~ /COSM/) && ($cov >= 5)) {
			print "Possibly somatic";
		} elsif (($af > 0.60) && ($af < 0.85) && ($cosm =~ /COSM/) && ($cov >= 5)) {
			print "Possibly somatic";
		} elsif (($af > 0.60) && ($af < 0.85) && ($rs =~ /rs/) && ($cosm !~ /COSM/) && ($cov >= 5)) {
			print "Possibly somatic";
		} elsif (($af < 0.4)  && ($ex_v eq '-') && ($cov >= 5)) {
			print "Somatic";
		} elsif (($af eq '.') || ($af eq '') || ($af eq '-')) {
			print "ND";
		} else {
			print "Unknown";
		}
		for ($y=24;$y<=$#a;$y++){
			$after = $after."\t".$a[$y];
		}
		print "$after";
	} else {
		$before = '';
		$after = '';
		@a = split/\t/;
		for ($i=0;$i<=$#a;$i++){
			$af_i = $i if ($a[$i] eq 'Sample AF');
			
			$dp_i = $i if ($a[$i] eq 'DP');
			
			$cosm_i = $i if ($a[$i] eq 'COSMIC_id');
			
			$rs_i = $i if ($a[$i] eq 'RS_id');
		}
		for ($y=0;$y<=23;$y++){
			$before = $before.$a[$y]."\t";
		}
		print "$before";
		print "Status";
		for ($y=24;$y<=$#a;$y++){
			$after = $after."\t".$a[$y];
		}
		print "$after";
	}
}
