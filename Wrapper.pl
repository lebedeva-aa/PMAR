use strict;
use warnings;
use Dir::Self;
use lib __DIR__ . '/../../lib'; 
use AtlasPipe;
use JSON;
use utf8;
use Getopt::Long 'GetOptions';
use Pod::Usage;
use Net::OpenSSH;
use Data::Dumper;
use Env;
use lib $ENV{AODADMIN};
use Aoddb;

#my %panel_DIC = (
#	"AODCPV1" => "CP",
#	"AODABCV1" => "ABC",
#	"CCP" => "CCP",
#	"CHPV2" => "CHP"
#	);

my $config		= __DIR__ . '/../../lib/Config.json';
my $panel_info		= __DIR__ . '/../../panel_info';
my $current_dir		= __DIR__;

my @tmpData;
my $tmpPath = (AtlasPipe::file_to_json($config))->{data_path}->{tmpPath};
my $samtools = (AtlasPipe::file_to_json($config))->{software}->{samtools};
my $hg19 = (AtlasPipe::file_to_json($config))->{data_path}->{hg19};
$config = AtlasPipe::file_to_json($config);


sub clearTmpData {
	foreach my $arg (@tmpData) {
		`rm -r $arg`;
		}
	@tmpData = ();
	}

sub fatality {
	clearTmpData;
	my $message = shift;
	print STDERR "$message\n" if defined $message;
	print STDERR "\nExit status 1\n";
	exit;
	}

sub indexBam {
	my $bam = shift;
	if (open(READ, "$bam.bai")) {
		close READ;
		} else {
		`$samtools index $bam`;
		}
	if (open(READ, "<$bam.bai")) {
		} else {die "Corrupted BAM file $bam\n"};
	}

sub run {
	my ($options) = @_;
	#die unless $options->{'mode'} eq 'parent';
	unless (defined($options->{'seed'})) {
		$options->{'seed'} = int(rand(100000000000000000));
		}
	print STDERR "Selected mode:\t",$options->{mode},"\n";
	print STDERR "Analysis seed:\t",$options->{seed},"\n";
	if ($options->{mode} eq 'parent') {
		$options->{remote} = remote_config($options);
		}
	my $bam = $options->{'b'};
#	my $panel = $panel_DIC{$options->{'panel'}};
	my $cellularity = $options->{'cell'};
	my $dir = $options->{'dir'};
	my $tmpPath = $config->{data_path}->{tmpPath} . "/";;
	my $folder = $tmpPath . "/" . $options->{seed};
	my $task = $current_dir . "/popa";;
	my $log_file = "$folder/run_log";
	if ($options->{mode} eq 'parent') {
		my $tmpPath = $options->{remote}->{data_path}->{tmpPath} . "/";
		my $bam_name = get_name($tmpPath . $options->{seed} . "." . get_base_name($bam));
		my $folder = $tmpPath . "/" . $options->{seed};
		my $task = $config->{remote}->{software_path} . "/Pipe/popa/popa";
		my $log_file = "$folder/run_log";
		create_remote_dir($options);
		} else {
		my $command = "mkdir -p -v $folder";
		print STDERR `$command`;
		}
	send_files($options, $bam, "$bam.bai");

	my $Analysis = $DB->Analysis($options->{a});
	my $dir = $Barcode->get_folder.$Analysis->get_id;
	`mkdir -v $dir`;
	my $command = "perl $task $bam_name $folder $panel $sex $cellularity > $log_file 2>&1";
	print STDERR "$command\n";
	my $response = get_response($command, $options);
	remove_tmp($options);
	system("perl $current_dir/post_analysis.pl ".$options->{dir}." ".$options->{b}. " $panel ". $options->{seed});
	}

sub worker {
        while ( my $passed = $work->dequeue ) {
                print STDERR "STARTED thread\n";
                my $path        = $passed->[0];
                my $name        = $passed->[1];
                my $fasta1      = $passed->[2];
                my $fasta2      = $passed->[3];

                `bwa mem -t 2 \$hg19 $fasta1 $fasta2 > $path/$name.sam`;
                `samtools view -@ 2 -b $path/$name.sam > $path/$name.bam`;
                `samtools sort -@ 2 -m 1G $path/$name.bam -o $path/$name.sorted.bam`;
                `mv $path/$name.sorted.bam $path/$name.bam`;
                `samtools index $path/$name.bam`;
                `rm $path/$name.sam`;

                }
        }
sub worker {
        while ( my $passed = $work->dequeue ) {
                my $bam         = $passed->[0];
                my $panel       = $passed->[1];
                my $path        = $passed->[2];
                my $name        = $passed->[3];
                `samtools mpileup -d 0 --reference \$hg19 $path/$name.bam > $path/$name.mpileup`;
                `perl /home/onco-admin/RnD/UEBAcall/mpileupToVcf.pl -input $path/$name.mpileup -output $path/$name.mpileup.vcf -limit 0`;
                `bgzip $path/$name.mpileup.vcf`;
                `tabix $path/$name.mpileup.vcf.gz`;
                `bcftools norm -cs -m -any -f \$hg19 $path/$name.mpileup.vcf.gz > $path/$name.mpileup.norm.vcf`;
                `mv $path/$name.mpileup.norm.vcf $path/$name.mpileup.vcf`;
                `rm $path/$name.mpileup $path/$name.mpileup.vcf.gz.tbi $path/$name.mpileup.vcf.gz`;

                `rm -r $path/analysis`;
                `perl /home/onco-admin/ATLAS_software/aod-pipe/Pipe/custom/recipe1.pl $bam $panel $path/`;
	}

sub main {
my $folder = $ARGV[0];
my $panel = $ARGV[1];
my $n_thread = $ARGV[2];
$n_thread = 10 unless defined $n_thread;

die "define path to folder with data with argument 1" unless defined $folder;
die "define panel Code with argument 2" unless defined $panel;
if (($panel ne "AODABCV1")and($panel ne "AODHRD15")and($panel ne "AODCPV2")and($panel ne "AODABCV2")) {
        die "unkown panel code";
        }

my $varList = __DIR__ . "/$ARGV[0]/Run.variant.list";
my $annotationVCF = __DIR__ . "/$ARGV[0]/Run.variant.annotation.vcf";
#my $annotationVCF = __DIR__ . "/$ARGV[0]/2.vcf";

`perl /home/onco-admin//RnD/UEBAcall/make_distributions.pl -l /home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/conf_data/${panel}_testkit_ILLMN/list_bam -v $annotationVCF -p /home/onco-admin/ATLAS_software/aod-pipe/panel_info/$panel/$panel.designed.bed -bdata /home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/conf_data/${panel}_testkit_ILLMN/bdata -mode append -n $n_thread`;

open (READ, "<$folder/runInfo");
open (WRITE, ">$folder/AODBETA.calls");

while (<READ>) {
        chomp;
        my @mas = split/\t/;
        next if $mas[0] =~ /#/;
        my $name = $mas[0];

        my $varList_local = "$folder/Fastq/$name/analysis/variant.list";
        my $vcf_input = "$folder/Fastq/$name/analysis/variant.vcf";
        `perl /home/onco-admin/ATLAS_software/aod-pipe/Prepare/varListToVCF.pl -l $varList_local -v $vcf_input`;

        #print STDERR "get counts\n";
        `perl /home/onco-admin//RnD/UEBAcall/get_counts_for_sample.pl -s $folder/Fastq/$name/$name.bam -v $vcf_input -o $folder/Fastq/$name/$name.cdata -p /home/onco-admin/ATLAS_software/aod-pipe/panel_info/$panel/$panel.designed.bed -n $n_thread`;
        #print STDERR "make_call\n";
        `perl /home/onco-admin//RnD/UEBAcall/make_call.pl -v $vcf_input -cdata $folder/Fastq/$name/$name.cdata -bdata /home/onco-admin/ATLAS_software/aod-pipe/Pipe/AODTKAbatch/conf_data/${panel}_testkit_ILLMN/bdata -n $n_thread > $folder/Fastq/$name/$name.calls`;

        #print STDERR "makeVCF\n";
        `perl /home/onco-admin//RnD/UEBAcall/makeVCF.pl -input $folder/Fastq/$name/$name.calls -output $folder/Fastq/$name/$name.AODBETA.vcf -sample $folder/Fastq/$name/$name.bam`;
        #print STDERR "filter VCF\n";
        `perl /home/onco-admin//RnD/UEBAcall/filterVCF.pl -input $folder/Fastq/$name/$name.AODBETA.vcf`;

        my $DB = AODDB->fast_connect();
        open (VCF, "<$folder/Fastq/$name/$name.AODBETA.vcf");

        while (<VCF>) {
                chomp;
                my $line = $_;
                next if m!#!;
                my @mas = split/\t/;
                my $mutation_name = "$mas[0]:$mas[1]$mas[3]>$mas[4]";
                my $Mutation = $DB->Mutation($mutation_name);
                my $rejectionStatus;
                my $r1 = $Mutation->isRejectedPFQ;
                my $r2 = $Mutation->isRejectedVCS;
                if ($r1 eq '0') {
                        if ($r2 eq '0') {
                                $rejectionStatus = 0;
                                } else {
                                $rejectionStatus = $r2;
                                }
                        } else {
                        $rejectionStatus = $r1;
                        }
                print WRITE "$name\t$rejectionStatus\t$mutation_name\t$line\n";
                }

        close VCF;
        }

close WRITE;
close READ;

}	


sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'	=> \$opts{'h'},
		'a|analysis=s'	=> \$opts{'b'}
	);
	return \%opts;
}

{
	my $options = option_builder();
	my $DB = AODDB->fast_connect();
	$options->{DB} = $DB;
	eval{run($options)};
	if ($@) {
		print STDERR "$@\n";
		die "Failed to run analyzer;\nExit status 1\n";
		}
}

__END__

=head1 NAME

=head1 SYNOPSIS

Main pipeline;

    -analysis  analysis code

=cut























