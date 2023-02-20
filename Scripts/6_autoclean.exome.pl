$in = $ARGV[0];
chomp $in;
open (IN, '<', $in) or die;
$in =~ s/\.txt/\.good\.txt/;
$good = $in;
$in =~ s/\.good\.txt/\.bad\.txt/;
$bad = $in;
open (GOOD, '>', $good) or die;
open (BAD, '>', $bad) or die;

while (<IN>) {
	if (/^Consequence/) {
		@a = split/\t/;
		for ($i=0;$i<=$#a;$i++) {
			$sample_f_i = $i if ($a[$i] eq 'Sample AF');
			$depth_i = $i if ($a[$i] eq 'DP');
			$mut_depth_i = $i if ($a[$i] eq 'Allele depth');
			$db_i = $i  if ($a[$i] eq 'DB');
			$rs_i = $i if ($a[$i] eq 'RS_id');
			$cosm_i = $i  if ($a[$i] eq 'COSMIC_id');
			$global_af_i = $i if ($a[$i] eq 'GLOBAL_AF');
			$eu_af_i = $i if ($a[$i] eq 'EU_AF');
			$global_gnomad_i = $i if ($a[$i] eq 'GLOBAL_AF_GNOMAD');
			$cons_i = $i if ($a[$i] eq 'Consequence');
			$name_i = $i if ($a[$i] eq 'Gene');
		}
		print GOOD "$_";
		print BAD "Reason of rejection\t$_";;
		next;
	}
	#print;
	@a = split/\t/;
	$mut_depth = $a[$mut_depth_i];
	$depth = $a[$depth_i];
	if (($a[$cosm_i] eq 'NA') && (($depth < 20) || ($a[$sample_f_i] < 0.05) || ($mut_depth < 10))) {
		print BAD "Nonconsistent with MSK recommendations for non-hotspots\t$_";
		#print "2\n";
	#non-hotspot
		next;
	}
	#previously found !!! ONLY FOR CCP
	if (($a[$db_i] > 0.9) && ($a[$rs_i] eq 'NA') && ($a[$cosm_i] eq 'NA')) {
		print BAD "Frequently occurred - possible artefact\t$_";
		#print "3\n";
		next;
	}
	#consequence
	if (($a[$cons_i] =~ /upstream_gene_variant/) || ($a[$cons_i] =~ /intron/) || ($a[$cons_i] =~ /synonymous/) || ($a[$cons_i] =~ /downstream_gene_variant/)) {	
		print BAD "Unwanted consequence\t$_";
		#print "5\n";
		next;
	}
	if (($a[$sample_f_i] >= 0.9) && ($a[$cons_i] =~ /frameshift_variant/)) {
		print BAD "Frameshift with high alleic frequency - possible sequencing error\t$_";
		#print "8\n";
		next;
	}
	print GOOD "$_";
}
