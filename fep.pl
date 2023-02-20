use strict;
use warnings;
use Env;
use lib $ENV{AODADMIN};
use Aoddb;
use Net::OpenSSH;
use Dir::Self;
use lib __DIR__ . '/../../lib';
use AtlasPipe;
use Data::Dumper;
use feature "switch";
no if $] >= 5.018, warnings => qw( experimental::smartmatch );

my $DB = AODDB->fast_connect();

# name of fastq files, directory

my $bam = $ARGV[0]; #enter two fastq.gz files separated with a comma
my $dir = $ARGV[1];
my $name = $ARGV[2];
$dir = "$dir/";
$name =~ s/\r?\n//g;
my $bam_marked = "$dir/$name.marked.bam";

my $current_dir	= __DIR__;
my $config	= "$current_dir/../../lib/Config.json";
$config		= AtlasPipe::file_to_json($config);

my %path;
$path{samtools}  = $config->{software}->{samtools};
$path{strelka}   = $config->{software}->{strelka};
$path{bcftools}  = $config->{software}->{bcftools};
$path{genome}    = $config->{data_path}->{genome};
$path{freebayes} = $config->{software}->{freebayes};
$path{picard}    = $config->{software}->{picard};
$path{bedtools}  = $config->{software}->{bedtools};
$path{pcgr}      = $config->{software}->{pcgr};
$path{pcgr_dir}  = $config->{software}->{pcgr_dir};
$path{pcgr_toml} = $config->{software}->{pcgr_toml};


my $strelka_dir = $dir.'Strelka';
Atlas::execute_cmd ("mkdir -p $strelka_dir");
my $pcgr_dir = $dir.'pcgr';
Atlas::execute_cmd ("mkdir -p $pcgr_dir");
#my $one_seq = `awk '{s++}END{print s/4}' $fastq1` if ($extension eq '.fastq');
#my $two_seq = `awk '{s++}END{print s/4}' $fastq2` if ($extension eq '.fastq');
#my $command = '| echo $((`wc -l`/4))';
#$one_seq = `zcat $fastq1 $command` if ($extension eq '.fastq.gz');
#$two_seq = `zcat $fastq2 $command` if ($extension eq '.fastq.gz');
#my $all_seq = $one_seq + $two_seq;
#system ("echo $all_seq > $size");
my $size = `$path{samtools} view -c $bam`;chomp $size;
my $sizeFile = $dir."/sequence.sizes.txt";
Atlas::execute_cmd ("echo $size > $sizeFile");

# Trimmomatic
#system ("java -jar /home/aod/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 4 $fastq1 $fastq2 $paired1 $unpaired1 $paired2 $unpaired2 ILLUMINACLIP:/home/aod/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50");
#system ("rm $unpaired1");
#system ("rm $unpaired2");

# samtools and Pickard
#system ("/home/aod/bin/samtools-1.9/samtools view -@ 4 -Sb $sam > $bam");
#system ("/home/aod/bin/samtools-1.9/samtools sort -@ 4 $bam > $bam_whole_sort");
#system ("/home/aod/bin/samtools-1.9/samtools index $bam_whole_sort");
#system ("/home/aod/bin/samtools-1.9/samtools view -@ 4 -h -L /home/aod/AOD/Clinical_exome/Scripts/agilent.onco.1.bed $bam > $sam_panel");
#system ("/home/aod/bin/samtools-1.9/samtools view -@ 4 -Sb $sam_panel > $bam_panel");
#system ("/home/aod/bin/samtools-1.9/samtools sort -@ 4 $bam_panel > $bam_sort");
my $bam_grep = "$dir/$name.grep.bam";
Atlas::execute_cmd ("$path{samtools} view -@ 4 -b -h -L $current_dir/conf/WES.bed $bam > $bam_grep");
my $metrics_file = "$dir/duplicate_metrics.txt";
#my $cmd = "java -jar $path{picard} MarkDuplicates I=$bam_grep O=$bam_marked M=$metrics_file REMOVE_DUPLICATES=true";
#print STDERR "$cmd\n";
#system ("$cmd");
Atlas::execute_cmd ("mv $bam_grep $bam_marked");
Atlas::execute_cmd ("$path{samtools} index $bam_marked");
#system ("rm $sam");
#system ("rm $bam");
#system ("rm $sam_panel");
#system ("rm $bam_panel");
#system ("rm $bam_sort");
#system ("rm $paired1");
#system ("rm $paired2");

# QC
my $depth = $dir.$name.'.depth.txt';
my $coverage = $dir.$name.'.coverage.txt';
my $table = $dir.$name.'.table.txt';
my $bad_genes = $dir.'/QC_data.txt';
my $poor_hotspots = $dir.'/poor_hotspots.txt';
my $uniformity = $dir.$name.'.uniformity.txt';
Atlas::execute_cmd ("$path{samtools} depth -d 0 -a -b $current_dir/conf/WES.bed $bam_marked > $depth");
Atlas::execute_cmd ("perl $current_dir/Scripts/get_coverage_new.pl $current_dir/conf/WES.bed < $depth > $coverage");
Atlas::execute_cmd ("perl $current_dir/Scripts/count_cov_exome.pl < $coverage > $table");
Atlas::execute_cmd ("perl $current_dir/Scripts/uniformity.pl $coverage < $depth > $uniformity");
Atlas::execute_cmd ("perl $current_dir/Scripts/get_bad_genes.pl $sizeFile $uniformity $table > $bad_genes");
Atlas::execute_cmd ("perl $current_dir/Scripts/4_poor_coverage_hotspots.new.pl $depth > $poor_hotspots");
#system ("rm $coverage");
#system ("rm $table");
#system ("rm $sizeFile");
#system ("rm $uniformity");

# variant calling (Strelka)
Atlas::execute_cmd ("$path{strelka} --bam $bam_marked --referenceFasta $path{genome} --runDir $strelka_dir");
my $starting_script = $strelka_dir.'/runWorkflow.py';
Atlas::execute_cmd ("$starting_script -m local");
my $vcf_gz = $strelka_dir.'/results/variants/variants.vcf.gz';
Atlas::execute_cmd ("gunzip $vcf_gz");
my $vcf = $strelka_dir.'/results/variants/variants.vcf';

# split  multiallelic Strelka VCF (MAXIM)
my $vcf_split = $strelka_dir.'/results/variants/variants.split.vcf';
Atlas::execute_cmd ("$path{bcftools} norm -m -any -f /home/aod/GENOME/hg19/hg19.fa $vcf > $vcf_split");
Atlas::execute_cmd ("mv $vcf_split $vcf");

# grep out not PASS variants
Atlas::execute_cmd ("grep -e \"\\sPASS\\s\" $vcf > $vcf.pass");
Atlas::execute_cmd ("grep '^#' $vcf > $vcf.header");
Atlas::execute_cmd ("cat $vcf.header $vcf.pass > $vcf");


# variant calling (Freebayes) (MAXIM)
Atlas::execute_cmd ("mkdir $dir/Freebayes");
Atlas::execute_cmd ("$path{freebayes} -f $path{genome} $bam_marked > $dir/Freebayes/result.vcf");
Atlas::execute_cmd ("$path{bcftools} norm -m -any -f $path{genome} $dir/Freebayes/result.vcf > $dir/Freebayes/result_split.vcf");
Atlas::execute_cmd ("perl $current_dir/Scripts/parse_FreeBayes.pl $dir/Freebayes/result_split.vcf > $dir/Freebayes/result_pass.vcf");
Atlas::execute_cmd ("grep '^#' $dir/Freebayes/result.vcf > $dir/Freebayes/header.vcf");
Atlas::execute_cmd ("sort -k1,1 -k2,2n -u $dir/Freebayes/result_pass.vcf | grep -v '^#' > $dir/Freebayes/result_uniq_content.vcf");
Atlas::execute_cmd ("cat $dir/Freebayes/header.vcf $dir/Freebayes/result_uniq_content.vcf > $dir/Freebayes/result_uniq.vcf");

# join VCF files (MAXIM)
# Ввиду вырожденности записи мутаций этот алгоритм может приводить к дублированной информации в итоговом VCF файле.
# В идеале нужен другой алгоритм конкатенации VCF файлов.
# В идеале необходим кастомный скрипт, который проверяет находятся ли две мутации в одном аллеле
Atlas::execute_cmd ("$path{bedtools} subtract -a $dir/Freebayes/result_uniq.vcf -b $vcf > $dir/Freebayes/result_addition.raw.vcf");
Atlas::execute_cmd ("sort -k1,1 -k2,2n -u $dir/Freebayes/result_addition.raw.vcf > $dir/Freebayes/result_addition.vcf");
Atlas::execute_cmd ("cat $vcf $dir/Freebayes/result_addition.vcf > $vcf.tmp");
Atlas::execute_cmd ("mv $vcf.tmp $vcf");

# change header
my $vcf_header = $dir.$name.'_header.vcf';
my $vcf_remake = "$dir/variant.filter.vcf";
Atlas::execute_cmd ("perl $current_dir/Scripts/fill_header.pl < $vcf > $vcf_header");
#Atlas::execute_cmd ("cp $vcf $vcf_header");
Atlas::execute_cmd ("perl $current_dir/Scripts/remake_vcf.pl < $vcf_header > $vcf_remake");
Atlas::execute_cmd ("rm $vcf_header");

Atlas::execute_cmd ("grep '^#' $vcf_remake > $vcf_remake.header");
Atlas::execute_cmd ("sort -k1,1 -k2,2n -u $vcf_remake | grep -v '^#' > $vcf_remake.content");
Atlas::execute_cmd ("cat $vcf_remake.header $vcf_remake.content > $vcf_remake");
Atlas::execute_cmd ("rm $vcf_remake.header $vcf_remake.content");


# pcgr
Atlas::execute_cmd ("$path{pcgr} --no_vcf_validate --input_vcf $vcf_remake $path{pcgr_dir} $pcgr_dir grch37 $path{pcgr_toml} $name");

# annotation
my $ann = $dir.'pcgr/'.$name.'.pcgr_acmg.grch37.snvs_indels.tiers.tsv';
my $vcf_gz_pcgr = $dir.'pcgr/'.$name.'.pcgr_acmg.grch37.pass.vcf.gz';
Atlas::execute_cmd ("gunzip $vcf_gz_pcgr");
my $vcf_pass = $dir.'pcgr/'.$name.'.pcgr_acmg.grch37.pass.vcf';
my $sort = $dir.$name.'.sort.txt';
my $gs = $dir.$name.'.sort.gs.txt';
my $marked = $dir.$name.'.sort.gs.marked.txt';
Atlas::execute_cmd ("perl $current_dir/Scripts/1_annotateCCP.new.pl $vcf_pass $ann Exome > $sort");
Atlas::execute_cmd ("perl $current_dir/Scripts/2_germ_som_est.pl < $sort > $gs");
Atlas::execute_cmd ("perl $current_dir/Scripts/3_mark_similar.pl < $gs > $marked");
Atlas::execute_cmd ("perl $current_dir/Scripts/6_autoclean.exome.pl $marked");
Atlas::execute_cmd ("rm $sort");
Atlas::execute_cmd ("rm $gs");

# cnv-calling
Atlas::execute_cmd ("perl $current_dir/Scripts/aod-cnv-we.pl $dir $bam $name");

