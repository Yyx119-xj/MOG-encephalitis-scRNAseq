#!/usr/bin/perl -w
my $ver="1.0.0";

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

#Before writing your programmeyou must write the detailed timediscriptionsparameter and it's explanation,Meanwhile,annotation your programme in English if possible.

my %opts;
GetOptions(\%opts,"gene=s","col=s","go=s", "sp=s", "kg=s","t=s", "pre=s", "c=s"  );

if(!defined($opts{gene}) || !defined($opts{go}) || !defined($opts{sp}) )
{
	print <<"	Usage End.";

	Description:This programme is used for

		Version: $ver

	Usage:perl $0

		-gene              infile                         must be given

		-col               int,column of gene             option, default is 1st

		-go                go des                         must be given

		-sp                species code(hsa,ath..) or "blast"              must be given

		-kg                kegg annot file                 must be given if "blast"

		-pre               prefix of outdir                option, default is "Gene"

		-t                 kegg anno type                  option, default is "gene:symbol"

		-c                 cutoff                          terms with less than cutoff number of genes are not used for statistical tests, default 5

	Usage End.

	exit;
}

###############Time_start##########
my $Time_Start;
$Time_Start = sub_format_datetime(localtime(time()));
print "\nStart Time :[$Time_Start]\n\n";
###################################

#Write your program
my $gene = $opts{gene};
my $column = defined($opts{col}) ? $opts{col} : 1;
$column = 1 if $column==0;
my $go_des = $opts{go};
my $species = $opts{sp};
if($species eq "blast" && !exists $opts{kg}) {
	die "No KEGG annotation file!\n";
}
my $kegg_des = defined $opts{kg} ? $opts{kg} : "";
my $prefix  = defined $opts{pre} ? $opts{pre} : "Gene";
# my $annotype  = defined $opts{t} ? $opts{t} : "gene:symbol";
my $annotype  = "gene:symbol";
my $cutoff = defined $opts{c} ? $opts{c} : 5;

my $current_dir = `pwd`;chomp($current_dir);

$gene = &AbsolutePath("file",$gene);
#$up_gene = &AbsolutePath("file",$up_gene);

&KEGG_pipe($gene,$species,$kegg_des,$prefix);
chdir($current_dir);
&GO_pipe($gene,$go_des,$prefix);

###############Time_end###########
my $Time_End;
$Time_End = sub_format_datetime(localtime(time()));
print "\nEnd Time :[$Time_End]\n\n";

# subroutines
sub KEGG_pipe {
	my ($g,$spe,$kg_des,$pre) = @_;
	my $outdir = $pre."_KEGG_analysis";
	`mkdir -p $outdir && cd $outdir`;
	chdir("$outdir");
	`cat $g | cut -f $column | sort | uniq | awk '\$1!="--" && \$1!="Gene" ' >sig_genelist`;
	#`cat $up | cut -f $column >up_genelist`;
	#`cat down_genelist up_genelist >sig_genelist`;
	#my $kg_des_cur = basename($kg_des);
	if ($spe eq "blast") {
		`ln -s -f $kg_des`;
		#`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/annotate.py -i $kg_des -t blastout:tab -s ko -o bg_annot`;
		`perl /public/pipeline/ablifebio/latest/GO_and_KEGG/KEGG/extract_fg_anno.pl -i $kg_des -g sig_genelist -o $pre\_KEGG_annotation.xls`;
		`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/identify.py -f $pre\_KEGG_annotation.xls -b $kg_des -d K -o $pre\_KEGG_pathway_iden.xls`;
	#	`perl /public/pipeline/ablifebio/latest/GO_and_KEGG/KEGG/generate_kegg_img.pl -up sig_genelist -down sig_genelist -iden $pre\_KEGG_pathway_iden.xls -anno $pre\_KEGG_annotation.xls -db /public/database/KEGG/kegg/pathway/ko -o $pre\_KEGG_pathway -imgdir img`;
	`sed -n '/^\\w/p' $pre\_KEGG_pathway_iden.xls > $pre\_KEGG_pathway.xls && sed -i '1i Term\tDatabase\tID\tInput number\tBackground number\tP-Value\tCorrected P-Value\tInput\tHyperlink' $pre\_KEGG_pathway.xls`;
		#`rm -rf $kg_des_cur`;
	} else {
		# print "python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/annotate.py -i sig_genelist -t $annotype -s $spe -o $pre\_KEGG_annotation.xls;python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/identify.py -f $pre\_KEGG_annotation.xls -b $spe -d K -o $pre\_KEGG_pathway_iden.xls;perl /public/pipeline/ablifebio/latest/GO_and_KEGG/KEGG/generate_kegg_img_genesymbol.pl -up sig_genelist -down sig_genelist -iden $pre\_KEGG_pathway_iden.xls -anno $pre\_KEGG_annotation.xls -species $spe -o $pre\_KEGG_pathway -imgdir img;\n";
		`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/annotate.py -i sig_genelist -t $annotype -s $spe -o $pre\_KEGG_annotation.xls`;
		# print "python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/annotate.py -i sig_genelist -t $annotype -s $spe -o $pre\_KEGG_annotation.xls\n";
		`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/identify.py -f $pre\_KEGG_annotation.xls -b $spe -d K -o $pre\_KEGG_pathway_iden.xls`;
	#	`perl /public/pipeline/ablifebio/latest/GO_and_KEGG/KEGG/generate_kegg_img_genesymbol.pl -up sig_genelist -down sig_genelist -iden $pre\_KEGG_pathway_iden.xls -anno $pre\_KEGG_annotation.xls -species $spe -o $pre\_KEGG_pathway -imgdir img`;
	`sed -n '/^\\w/p' $pre\_KEGG_pathway_iden.xls > $pre\_KEGG_pathway.xls && sed -i '1i Term\tDatabase\tID\tInput number\tBackground number\tP-Value\tCorrected P-Value\tInput\tHyperlink' $pre\_KEGG_pathway.xls`;
	}
	# top term bar plot
	#`cut -f 1-7 $pre\_KEGG_pathway.xls > kegg.txt && Rscript /public/pipeline/ablifegraphic/latest/ablifeR/pipeline/bar_gokegg_term.r -f kegg.txt -n $pre\_KEGG_top`;
	`cut -f 1-7 $pre\_KEGG_pathway.xls > kegg.txt && Rscript /public/pipeline/ablifepipe/latest/src/SingleCell/pipeline_v3/Scatter_multiSample.r -f kegg.txt -n $pre\_KEGG_top_Scatter -t $pre\_KEGG_top_Scatter`;
	#`head -n 11 $pre\_KEGG_pathway.xls | cut -f 1-7 > $pre\_KEGG_pathway && Rscript /public/pipeline/ablifegraphic/latest/ablifeR/pipeline/gokegg_Scatter_multiSample.r -f $pre\_KEGG_pathway -n $pre\_KEGG_top -t $pre\_KEGG_top10`;

}

sub GO_pipe {
	my ($g,$go,$pre) = @_;
	my $outdir = $pre."_GO_analysis";
	`mkdir -p $outdir && cd $outdir`;
	chdir("$outdir");
	`ln -s -f $go`;
	my $go_cur = basename($go);
	`cat $g | cut -f $column | sort | uniq | awk '\$1!="--" && \$1!="Gene" ' >sig_genelist`;
	#`cat $up | cut -f $column >up_genelist`;
	#`cat down_genelist up_genelist >sig_genelist`;
	`perl /public/pipeline/ablifebio/latest/GO_and_KEGG/GO/extract_fg_anno_go.pl -i $go_cur -g sig_genelist -o $pre\_GO_annotation`;
	`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/goanalyze.py -f $pre\_GO_annotation_C.xls -b $go_cur\_C -o $pre\_GO_enrichment_C.xls -c $cutoff && perl /public/pipeline/ablifebio/latest/GO_and_KEGG/GO/generate_go_web.pl -up sig_genelist -down sig_genelist -iden $pre\_GO_enrichment_C.xls -anno $pre\_GO_annotation_C.xls -o $pre\_GO_enrichment_web_C`;
	`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/goanalyze.py -f $pre\_GO_annotation_F.xls -b $go_cur\_F -o $pre\_GO_enrichment_F.xls -c $cutoff && perl /public/pipeline/ablifebio/latest/GO_and_KEGG/GO/generate_go_web.pl -up sig_genelist -down sig_genelist -iden $pre\_GO_enrichment_F.xls -anno $pre\_GO_annotation_F.xls -o $pre\_GO_enrichment_web_F`;
	`python2.7 /public/software/biotools/KOBAS/kobas2.0-20131201/scripts/goanalyze.py -f $pre\_GO_annotation_P.xls -b $go_cur\_P -o $pre\_GO_enrichment_P.xls -c $cutoff && perl /public/pipeline/ablifebio/latest/GO_and_KEGG/GO/generate_go_web.pl -up sig_genelist -down sig_genelist -iden $pre\_GO_enrichment_P.xls -anno $pre\_GO_annotation_P.xls -o $pre\_GO_enrichment_web_P`;
	`rm -rf $go_cur* `;
	# top term bar plot
	#`cat $pre\_GO_enrichment_C.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_C_top.txt && Rscript /public/pipeline/ablifegraphic/latest/ablifeR/pipeline/bar_gokegg_term.r -f $pre\_GO_enrichment_C_top.txt -n $pre\_GO_enrichment_C_top`;
	#`cat $pre\_GO_enrichment_P.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_P_top.txt && Rscript /public/pipeline/ablifegraphic/latest/ablifeR/pipeline/bar_gokegg_term.r -f $pre\_GO_enrichment_P_top.txt -n $pre\_GO_enrichment_P_top`;
	#`cat $pre\_GO_enrichment_F.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_F_top.txt && Rscript /public/pipeline/ablifegraphic/latest/ablifeR/pipeline/bar_gokegg_term.r -f $pre\_GO_enrichment_F_top.txt -n $pre\_GO_enrichment_F_top`;
	
	`cat $pre\_GO_enrichment_C.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_C_top.txt && Rscript /public/pipeline/ablifepipe/latest/src/SingleCell/pipeline_v3/Scatter_multiSample.r -f $pre\_GO_enrichment_C_top.txt -n $pre\_GO_enrichment_C_top_Scatter -t $pre\_GO_enrichment_C_top_Scatter`;
	`cat $pre\_GO_enrichment_P.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_P_top.txt && Rscript /public/pipeline/ablifepipe/latest/src/SingleCell/pipeline_v3/Scatter_multiSample.r -f $pre\_GO_enrichment_P_top.txt -n $pre\_GO_enrichment_P_top_Scatter -t $pre\_GO_enrichment_P_top_Scatter`;
	`cat $pre\_GO_enrichment_F.xls | perl -ne 'chomp;if(/^\\#Term/){\$i++;}if(\$i==3){\$flag=1;}if(\$flag==1){next if(/^\\s*\$/);next if(/^----/);\$_=~s/^\\#//;print \$_,"\\n";}' > $pre\_GO_enrichment_F_top.txt && Rscript /public/pipeline/ablifepipe/latest/src/SingleCell/pipeline_v3/Scatter_multiSample.r -f $pre\_GO_enrichment_F_top.txt -n $pre\_GO_enrichment_F_top_Scatter -t $pre\_GO_enrichment_F_top_Scatter`;
}

sub AbsolutePath {		# Get the absolute path of the target directory or file
	my ($type,$input) = @_;
	my $return;
	if ($type eq "dir"){
		my $pwd = `pwd`;
		chomp $pwd;
		chdir($input);
		$return = `pwd`;
		chomp $return;
		chdir($pwd);
	} elsif($type eq 'file') {
		my $pwd = `pwd`;
		chomp $pwd;
		my $dir=dirname($input);
		my $file=basename($input);
		chdir($dir);
		$return = `pwd`;
		chomp $return;
		$return .="\/".$file;
		chdir($pwd);
	}
	return $return;
}

###############Sub_format_datetime
sub sub_format_datetime {#Time calculation subroutine
    my($sec, $min, $hour, $day, $mon, $year, $wday, $yday, $isdst) = @_;
	$wday = $yday = $isdst = 0;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}
