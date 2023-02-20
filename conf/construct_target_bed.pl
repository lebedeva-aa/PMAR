use strict;
use warnings;
use Env;
use lib $ENV{AODADMIN};
use Aoddb;
use Net::OpenSSH;
use Dir::Self;
use lib __DIR__ . '/../../../lib';
use AtlasPipe;
use Data::Dumper;
use feature "switch";
no if $] >= 5.018, warnings => qw( experimental::smartmatch );

my $current_dir = __DIR__;
my $config      = "$current_dir/../../../lib/Config.json";
$config         = AtlasPipe::file_to_json($config);
my %path;
$path{GTF} = $config->{data_path}->{GTF};
$path{clinvar} = $config->{data_path}->{clinvar};


`rm WES.bed`;
`rm WES.clinvar.bed`;
open (GENE_LIST, "<gene_list");

while (<GENE_LIST>) {
	chomp;
	my $gene_symbol = $_;
	my $grep_cmd = "grep \"gene_name \\\"$gene_symbol\\\"\" $path{GTF} | grep \"\\stranscript\\s\" | grep \"tag \\\"basic\\\"\" | grep \"transcript_biotype \\\"protein_coding\\\"\" | grep \"tag \\\"CCDS\\\"\"";
	my $transcript_output = `$grep_cmd`;chomp $transcript_output;
	my @transcript = split/\n/, $transcript_output;
	print STDERR "$gene_symbol\n";
	if ((scalar @transcript) eq 0) {
		$grep_cmd = "grep \"gene_name \\\"$gene_symbol\\\"\" $path{GTF} | grep \"\\stranscript\\s\" | grep \"tag \\\"basic\\\"\" | grep \"transcript_biotype \\\"protein_coding\\\"\"";
		$transcript_output = `$grep_cmd`;chomp $transcript_output;
		@transcript = split/\n/, $transcript_output;
		}
	die "No transcripts found for gene $gene_symbol\n" if (scalar @transcript) eq 0;
	foreach my $arg (@transcript) {
		my $transcript_id = grep_field($arg, "transcript_id");
		if ($transcript_id) {
			print STDERR "$gene_symbol\t$transcript_id\n";
			$grep_cmd = "grep \"transcript_id \\\"$transcript_id\\\"\" $path{GTF} | grep \"\\sCDS\\s\"";
			#my $line = `$grep_cmd`;chomp $line;
			foreach my $line (split/\n/, `$grep_cmd`) {
				my $exon = grep_field($line, "exon_number");
				`echo '$line' | awk '{min = (\$4 < \$5) ? \$4 : \$5; max = (\$4 < \$5) ? \$5 : \$4; print \"chr\"\$1 \"\\t\" min-5 \"\\t\" max+5 \"\\t$gene_symbol:$transcript_id:$exon\"}' >> WES.bed`;
				}
			} else {
			die;
			}
		}
	my $cmd = "grep -E \"GENEINFO=[^;]+\\|$gene_symbol:|GENEINFO=$gene_symbol:\" $path{clinvar} | grep \"Pathogenic\\|Likely_pathogenic\" | awk '{print \"chr\"\$1 \"\\t\" \$2-1 \"\\t\" \$2+1 \"\\tCV_\" \$3}' >> WES.clinvar.bed";
	`$cmd`;
	next;
	}

close GENE_LIST;
`sort -k1,1 -k2,2n -u WES.bed > WES.sorted.bed`;
`bedtools merge -i WES.sorted.bed -d 5 -c 4 -o collapse > WES.bed`;
`cp WES.bed WES.exome.bed`;

`bedtools subtract -a WES.clinvar.bed -b WES.bed > WES.clinvar.addition.bed`;
`cat WES.bed WES.clinvar.addition.bed > WES.clinical.bed`;
`mv WES.clinical.bed WES.bed`;

`sort -k1,1 -k2,2n -u WES.bed > WES.sorted.bed`;
`bedtools merge -i WES.sorted.bed -d 5 -c 4 -o collapse > WES.bed`;

`sed -e 's/:.*//' WES.exome.bed > WES.exome.gene.bed`;

`rm WES.sorted.bed WES.clinvar.addition.bed WES.clinvar.bed`;

sub grep_field {
	my $line = shift;
	my $field = shift;
	if ($line =~ /$field "([^"]+)"(.*)/) {
		return $1;
		} else {
		return undef;
		}
	}


