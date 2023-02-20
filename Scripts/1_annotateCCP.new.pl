use Dir::Self;
use Data::Dumper;
use lib __DIR__ . '/../../lib';
my $current_dir = __DIR__;

use Env;
use lib $ENV{AODADMIN};
use Aoddb;

my $DB = AODDB->fast_connect();

open (READ1, "<$ARGV[0]"); #input VCF from Ion Torrent
open (READ2, "<$ARGV[1]"); #annotation from VEP
$panel_type = $ARGV[2];
open (EXONS, '<', "$current_dir/hg19.exons") or die;

%exons = ();
while (<EXONS>) {
	s/\r?\n//g;
	@sp = split/\t/;
	$exons{$sp[0]} = $sp[1];
}
close EXONS;

#%freq;
#%depth;
%countRef = ();
#%count;
$conseq = '';
$transcr = '';
$hgvsc = '';
$exon = '';
$hgvsp = '';
$number = 0;

#open (LEGEND, "<CCPlegend.info");
$dir = '/home/aod/ATLAS/CCP/' if ($panel_type eq 'CCP');
$dir = '/home/aod/ATLAS/ABC/' if ($panel_type eq 'ABC');
$dir = '/home/aod/ATLAS/CHP/' if ($panel_type eq 'CHP');
$dir = '/home/aod/ATLAS/CCP/' if ($panel_type eq 'CP');
$dir = '/home/aod/ATLAS/RCMGLYNCHV1/' if ($panel_type eq 'RCMGLYNCHV1');
$dir = '/home/aod/ATLAS/Exome/' if ($panel_type eq 'Exome');

opendir DIR, $dir or die;

while ($file = readdir DIR) {
	next unless ($file =~ /\.vcf$/);
	$number += 1;
	$file = $dir.$file;
	open (DATA, '<', $file) or die;
	while (<DATA>) {
		next if m!^#!;
		@data = split/\t/;
		@dalt = split/,/, $data[4];
		foreach $darg (@dalt) {
			++$countRef{"$data[0]:$data[1]$data[3]>$darg"};
			}
		}
	}
#close LEGEND;
#@vcf_file = <READ1>;
@annotation = <READ2>;
%vcf_file = ();
%global = ();
%eu = ();
%global_gnomad = ();
$number = 1 if $number eq 0;
$number = "SELECT COUNT(*) FROM Analysis INNER JOIN Barcode ON Barcode.barcodeName = Analysis.barcodeName WHERE panelCode = 'WES' and analysisRole = 'Major';";
$number = $DB->execute_select_single($number);
print STDERR "NUMBER OF SAMPLES: $number\n";
while (<READ1>) {
	next if (/^#/);
	$af = '.';
	$dp = '.';
	$allele_depth = '.';
	$saf = '.';
	$sar = '.';
	$read_distr = '.';
	$count = '.';
	@mas = split/\t/;
	my $var_name = lc("$mas[0]:$mas[1]$mas[3]>$mas[4]");
	$var_name = "chr$var_name" unless $var_name =~ /^chr/;
	$ref = $mas[3];
	my $chrom = $mas[0];
	$chrom = "chr$chrom" unless $chrom =~ /^chr/;
	#$count = $countRef{'chr'."$mas[0]:$mas[1]$mas[3]>$mas[4]"} if (exists($countRef{'chr'."$mas[0]:$mas[1]$mas[3]>$mas[4]"}));
	my $sql = "SELECT COUNT(*) FROM MutationResult INNER JOIN Mutation ON Mutation.mutationId = MutationResult.mutationId INNER JOIN Analysis ON Analysis.analysisName = MutationResult.analysisName INNER JOIN Barcode ON Barcode.barcodeName = Analysis.barcodeName WHERE mutationGenomicPos = '$mas[1]' AND mutationChr = '$chrom' AND mutationRef = '$mas[3]' AND mutationAlt = '$mas[4]' AND panelCode = 'WES' and analysisRole = 'Major';";
	print STDERR "$sql\n";
	$count = $DB->execute_select_single($sql);
	$count = $count/$number;
	$line = $mas[0].':g.'."$mas[1]$mas[3]>$mas[4]";
	#$line =~ s/^chr//;
	$global{$line} = 'NA';
	$eu{$line} = 'NA';
	$global_gnomad{$line} = 'NA';
	@info = split/;/, $mas[7];
	foreach $q (@info) {
		if ($q =~ /^AF=/) {
			$af = $q;
			$af =~ s/^AF=//;
		} elsif ($q =~ /^DP=/) {
			$dp = $q;
			$dp =~ s/^DP=//;
		} elsif ($q =~ /^AO/) {
			$allele_depth = $q;
			$allele_depth =~ s/^AO=//;
		} elsif ($q =~ /^SAF=/) {
			$saf = $q;
			$saf =~ s/^SAF=//;
		} elsif ($q =~ /^SAR=/) {
			$sar = $q;
			$sar =~ s/^SAR=//;
		} elsif ($q =~ /^GLOBAL_AF_1KG=/) {
			$global{$var_name} = $q;
			$global{$var_name} =~ s/^GLOBAL_AF_1KG=//;
		} elsif ($q =~ /^EUR_AF_1KG=/) {
			$eu{$var_name} = $q;
			$eu{$var_name} =~ s/^EUR_AF_1KG=//;
		} elsif ($q =~ /^GLOBAL_AF_GNOMAD=/) {
			$global_gnomad{$var_name} = $q;
			$global_gnomad{$var_name} =~ s/^GLOBAL_AF_GNOMAD=//;
		}
	}
	$read_distr = $saf.' : '.$sar;
	$vcf_file{$var_name} = "$af\t$dp\t$allele_depth\t$read_distr\t$count\t";
}
print STDERR Dumper \%vcf_file;
my %collumnDic;
my $i = 0;
map {$collumnDic{$_} = $i; ++$i} (split/\t/,$annotation[0]);
print "Consequence\tGene\tExon\tChromosome\tDP\tSample AF\tHGVSc\tHGVSp\tDB\tAllele depth\tRead distribution\tStrand bias\tRS_id\tCOSMIC_id\tGLOBAL_AF\tEU_AF\tGLOBAL_AF_GNOMAD\tTIER\tCLINVAR\tCLINVAR_CLNSIG\tONCOGENE\tTUMOR_SUPPRESSOR\tEFFECT_PREDICTIONS\tCANCER_MUTATION_HOTSPOT\tINTOGEN_DRIVER_MUT\tTCGA_PANCANCER_COUNT\tTCGA_FREQUENCY\tICGC_PCAWG_OCCURRENCE\tCHEMBL_COMPOUND_ID\tCHEMBL_COMPOUND_TERMS\tENSEMBL_TRANSCRIPT_ID\tPROTEIN_DOMAIN\t$annotation[0]";
for ($i=1;$i<=$#annotation;$i++) {
	@mas2 = split/\t/, $annotation[$i];
	$chrom = "$mas2[0]:$mas2[1]$mas2[2]>$mas2[3]";
	$chrom = 'Chr'.$chrom;
	my $var_name = lc($chrom);
	$chrom =~ s/g\.//;
	$gene = $mas2[8];
	if ($mas2[$collumnDic{"CDS_CHANGE"}] eq 'NA') {
		$exon = '-';
		if ($mas2[23] ne 'NA') {
			$hgvsc = $mas2[23];
			@s1 = split/:/, $hgvsc;
			$hgvsc = $s1[1];
		} else {
			$hgvsc = 'NA';
		}
		$hgvsp = '-';
	} else {
		($conseq, $transcr, $hgvsc, $exon, $hgvsp) = split/:/, $mas2[$collumnDic{"CDS_CHANGE"}];
		@q = split/\./, $transcr;
		$transcr = $q[0];
		if (($exon eq '') || ($exon eq 'NA')) {
			$exon = '-';
		} else {
			$exon = $exon.'/'.$exons{$transcr};
			$exon =~ s/^exon//;
		}
	}
	$hgvsp = '-' if ($hgvsp eq '');
	$hgvsc = '-' if ($hgvsc eq '');
	$conseq = $mas2[$collumnDic{"CONSEQUENCE"}];
	$rs_id = $mas2[$collumnDic{"DBSNPRSID"}];
	$cosmic_id = $mas2[$collumnDic{"COSMIC_MUTATION_ID"}];
	$cosmic_id =~ s/&/,/g;
	$clinvar = $mas2[$collumnDic{"CLINVAR"}];
	$clinvar_clnsig = $mas2[$collumnDic{"CLINVAR_CLNSIG"}];
	$protein_domain = $mas2[$collumnDic{"PROTEIN_DOMAIN"}];
	$oncogene = $mas2[$collumnDic{"ONCOGENE"}];
	$tsupressor = $mas2[$collumnDic{"TUMOR_SUPPRESSOR"}];
	$tier = $mas2[$collumnDic{"TIER"}];
	$eff_pred = $mas2[$collumnDic{"EFFECT_PREDICTIONS"}];
	$cmh = $mas2[$collumnDic{"MUTATION_HOTSPOT"}];
	$idm = "-";
	@spkl = split/\t/, $vcf_file{$var_name};
	print STDERR $var_name,"\n";
	@bo = split/\s/, $spkl[3];
	$strand = '';	
	if (($bo[0] == 0) || ($bo[2] == 0)) {
		$strand = 'BAD';
	} elsif (($bo[0] < $bo[2]) && ($bo[2]/$bo[0] < 2)) {
		$strand = 'GOOD';
	} elsif (($bo[2] < $bo[0]) && ($bo[0]/$bo[2] < 2)) {
		$strand = 'GOOD';
	} elsif ($bo[2] == $bo[0]) {
		$strand = 'GOOD';
	} else {
		$strand = 'BAD';
	}
#	print STDERR "$vcf_file{$mas2[0]}\n";
	#print "$chrom\t$gene\t$conseq\t$transcr\t$exon\t$hgvsc\t$hgvsp\t$rs_id\t$cosmic_id\t$clinvar\t$protein_domain\t$oncogene\t$tsupressor\t$tier\t$global{$mas2[0]}\t$eu{$mas2[0]}\t$annotation[$i]";
	print "$conseq\t$gene\t$exon\t$chrom\t$spkl[1]\t$spkl[0]\t$hgvsc\t$hgvsp\t$spkl[4]\t$spkl[2]\t$spkl[3]\t$strand\t$rs_id\t$cosmic_id\t$global{$var_name}\t$eu{$var_name}\t$global_gnomad{$var_name}\t$tier\t$clinvar\t$clinvar_clnsig\t$oncogene\t$tsupressor\t$eff_pred\t$cmh\t$idm\t$transcr\t$protein_domain\n";
}

close READ1;
close READ2;


