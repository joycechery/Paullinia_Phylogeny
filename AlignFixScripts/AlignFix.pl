########################################################################################################
########################################################################################################
# script to align Dimocarpus lognan transcriptome reads to Citrus Genome, Call SNPs and make consensus sequences closer to Sapindaceae#
#  written by Chodon Sass chodon at gmail.com, 2015
########################################################################################################

#!/Usr-/bin/perl -w
 use warnings;
use strict;
use Bio::SeqIO;


##define libraries, home directory and extension to home directory plus library where bams and vcfs are

my @lib = qw(JCSR);

my $novoalign = '/global/home/users/chodon/bin/novocraft/novoalign';
my $novoindex = '/global/home/users/chodon/bin/novocraft/novoindex';
my $gatk = '/global/home/users/chodon/bin/GenomeAnalysisTK.jar';
my $markDups = '/global/home/users/chodon/bin/MarkDuplicates.jar';
my $createDict = '/global/home/users/chodon/bin/CreateSequenceDictionary.jar';
my $readGroups = '/global/home/users/chodon/bin/AddOrReplaceReadGroups.jar';

my $home = '/clusterfs/rosalind/users/chodon/joyce/';
my $selectseqs = 'selectseqs.pl'; 

###### start going through libraries

foreach my $lib (@lib){
    my $dir = $home;
    my $tempDir = $dir . 'temp/';
    mkdir ($tempDir) unless -d ($tempDir);
    my $clean1 = $dir . $lib . '_1_final.txt';
    my $clean2 = $dir . $lib . '_2_final.txt';
    my $cleanu = $dir . $lib . '_u_final.txt';
    my $refs = $dir . 'refs_JGC.fa';
     my $dict = substr ($refs, 0, -2) . 'dict';
    my $indexed_assemblies_in_target =  substr ($refs, 0, -2) . "nix";
    my $bam =  $dir . $lib . '_m2n240_4.bam';
  system("samtools faidx $refs");
  system("$novoindex $indexed_assemblies_in_target $refs");
  system("$novoalign -d $indexed_assemblies_in_target -f $clean1 $clean2 -i PE 230, 70 -t 480 -r Random -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outPairedSam1");
  system("$novoalign -d $indexed_assemblies_in_target -f $cleanu -t 480 -r Random -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outSoloSam1");
  system("grep -v \'ZS:Z:NM\' outPairedSam1 > target_pair.sam");
  system("grep -v \'ZS:Z:NM\' outSoloSam1 > target_solo.sam");
  system("samtools view -bS target_pair.sam > target_pair.bam");
  system("samtools view -bS target_solo.sam > target_solo.bam");
  system("samtools merge -f target.bam target_solo.bam target_pair.bam");
  system("samtools sort target.bam target.sorted");
  system("samtools index target.sorted.bam");
  system("java -Djava.io.tmpdir=$tempDir -jar $markDups INPUT=target.sorted.bam OUTPUT=target.duped.bam METRICS_FILE=target.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
  system("samtools sort target.duped.bam target.newduped");
  system("samtools index target.newduped.bam");
  system("java -Djava.io.tmpdir=$tempDir -jar $createDict REFERENCE=$refs OUTPUT=$dict");
  system("java -jar $readGroups INPUT=target.newduped.bam OUTPUT=$bam RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
  system("samtools index $bam");
system("rm target.duped.bam target.sorted.bam target.genome target.sorted.bam.bai target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.metric target.bam target.newduped.bam target.newduped.bam.bai");

    my $newcov = $dir . $lib . '.newcov';
    my $vcf = $dir . $lib . '.filter.vcf';
    my $filterQual = q("QD < 2.0 || HaplotypeScore > 13.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0");
    my $filterDp20 = '"DP < 20 && DP > 10"';
    my $filterDp10 = '"DP < 11"';
    my $filterHethigh = '"ABHet > 0.699"';
    my $filterHetlow = '"ABHet < 0.200"';
    my $filterAltpos = '"ABHet < 0.500 && ABHet > 0.199"';
    my $filterAltdef = '"ABHom == 1.00"';
    my $filterHomhigh = '"ABHom > 0.699 && ABHom < 1.00"';
    my $filterHomlow = '"ABHom > 0.499 && ABHom < 0.700"';
    my $exp = '--filterExpression';
    my $name = '--filterName';
    my $intervals = substr ($refs, 0, -2) . "intervals";
   system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T RealignerTargetCreator -R $refs -I $bam -o $intervals");
   system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T IndelRealigner -R $refs -I $bam -targetIntervals $intervals -o $bam.realigned.bam -LOD 0.1 -model USE_SW");
   system("samtools index $bam.realigned.bam");
   system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T HaplotypeCaller -R $refs -I $bam.realigned.bam -dt NONE -o $dir$lib.hc.vcf");
   system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T ReadBackedPhasing -R $refs -I $bam.realigned.bam --variant $dir$lib.hc.vcf --min_base_quality_score 10 -dt NONE -o $dir$lib.raw.vcf");
   system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T VariantAnnotator -A DepthPerAlleleBySample -A AlleleBalance -A FisherStrand -A HaplotypeScore -A HardyWeinberg -R $refs -I $bam.realigned.bam --variant $dir$lib.raw.vcf -o $dir$lib.annotated.vcf");
   system("java -Djava.io.tmpdir=$tempDir -jar $gatk -T VariantFiltration -R $refs -o $vcf --variant $dir$lib.annotated.vcf $exp $filterQual $name \"qual\" $exp $filterDp20 $name \"dp20\" $exp $filterDp10 $name \"dp10\" $exp $filterHethigh $name \"hethigh\" $exp $filterHetlow $name \"hetlow\" $exp $filterAltpos $name \"altpos\" $exp $filterAltdef $name \"altdef\" $exp $filterHomhigh $name \"homhigh\" $exp $filterHomlow $name \"homlow\"");

#name vcf file
   
    my $rawvcf = $dir. $lib . '.filter.vcf';
    my $processvcf2 = $dir . $lib . '.short1.vcf';
    my $columnName = $dir . $lib . 'columns';
#this annotates indels and when there is more than one base as an alternative
    system("grep -v \'\#\' $rawvcf | awk \'\$5 !~ /[A-Z],[A-Z]/\'  \| sed \'/ABH/!s/PASS/INDEL/\' \| sed \'/ABH/!s/dp10/LOWDEL/\' \| sed \'/ABH/!s/dp20/LOWDEL/\' > $processvcf2");
    system("grep \'\#CHROM\' $rawvcf > $columnName");

    my $vcffil = $dir. $lib . '.short1.vcf';
    my $vcfv2 = $dir. $lib . '.short2.vcf';
    my $vcfv3 = $dir. $lib . '.short3.vcf';

### this goes through gatk vcf and filters variants for printing onto a second vcf that will be used to generate a new fasta for gatk again -- can't have variants
   open(VCF, "<$vcffil");
    open(NEWVCF, ">>$vcfv2");

foreach my $vcline (<VCF>)    {
       
       my @tmp = split(/\t/,$vcline);
       my $prea0; my $prea1; my $join; my $newletter; my $newline;
       my $a0; my $a1;

### for coverage less than 10, use 50% rule unless quality, then leave ref
       if ($tmp[6] =~ /dp10/){
          if ($tmp[6] =~ m/altdef|altpos|homhigh|homlow|hetlow/) {
               if ($tmp[6] !~ m/qual/) {
                  print NEWVCF $vcline;
               }
           }                                                                                                                                                  
      }


# for coverage less than 20, use 50% rule, unless qual, then leave ref
             
if ($tmp[6] =~ m/dp20/) {
   if ($tmp[6] =~  m/altdef|hetlow|homhigh|homlow|altpos/) {
print NEWVCF $vcline;
    } #print
             
}
# for qual, use 50 percent rule (this is vestigial, could be incorporated
               
if ($tmp[6] =~ m/qual/) {
 if ($tmp[6] !~ m/dp20|dp10/) {
if ($tmp[6] =~ m/hetlow|altpos|altdef|homhigh|homlow/) {
   print NEWVCF $vcline;
}# print        
   }              
}
## no dp and no qual

if ($tmp[6] !~ m/qual|dp20|dp10/) {
    if ($tmp[6] =~ m/hetlow|homhigh|altdef|altpos|homlow/) {
print NEWVCF $vcline;
    } #print
} 

## deal with indels, if INDEL is >20x coverage and 50% rule, print to new ref

if ($tmp[6] =~ m/INDEL/) {
    my @tmpindel = split(/:/,$tmp[9]);
    my $indelcov = $tmpindel[1];
    my @tmpindelcov = split(/,/,$indelcov);
    if (($tmpindelcov[1] + $tmpindelcov[0]) > 20) {
if ($tmpindelcov[1] > $tmpindelcov[0]) {
    print NEWVCF $vcline;
}
    }
}
             
###
    } 
   my $newFasta =  $dir. $lib . '.newFasta.fa';
   my $getnoemptylines = "\'\^\$\'";
    system("cat $columnName $vcfv2 | grep -v $getnoemptylines > $vcfv3");
   my $zipped = $vcfv3 . '.gz';
    system("bgzip $vcfv3");
    system("tabix -p vcf $zipped");
    system("cat $refs | vcf-consensus $zipped > $newFasta");
print "now finished with vcf - going to run second novoalign";


# run novoalign
 
   
    my $insertSize = 240;
    my $insertStd = 70;
    my $maxScore = 90;
    my $nref = $dir .  $lib . '.newFasta.fa';
    my $out1 =  $dir .  $lib . '.new2';
    my $out2 =  $dir . $lib . '.new2.bam';
    my $indexed_assemblies_in_target2 =  substr ($nref, 0, -2) . "nix";
    my $dict2 = substr ($nref, 0, -2) . 'dict';
  system("samtools faidx $nref");
  system("$novoindex $indexed_assemblies_in_target2 $nref");
  system("$novoalign -d $indexed_assemblies_in_target2 -f $clean1 $clean2 -i PE $insertSize, $insertStd -t $maxScore -r Random -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outPairedSam2");
  system("$novoalign -d $indexed_assemblies_in_target2 -f $cleanu -t $maxScore -r Random -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -F STDFQ -o SAM > outSoloSam2");
  system("grep -v \'ZS:Z:NM\' outPairedSam2 > target2_pair.sam");
  system("grep -v \'ZS:Z:NM\' outSoloSam2 > target2_solo.sam");
  system("samtools view -bS target2_pair.sam > target2_pair.bam");
  system("samtools view -bS target2_solo.sam > target2_solo.bam");
  system("samtools merge -f target2.bam target2_solo.bam target2_pair.bam");
  system("samtools sort target2.bam target2.sorted");
  system("samtools index target2.sorted.bam");
  system("java -Djava.io.tmpdir=$tempDir -jar $markDups INPUT=target2.sorted.bam OUTPUT=target2.duped.bam METRICS_FILE=target2.metric REMOVE_DUPLICATES=true ASSUME_SORTED=true");
  system("samtools sort target2.duped.bam target2.newduped");
  system("samtools index target2.newduped.bam");
  system("java -Djava.io.tmpdir=$tempDir -jar $createDict  REFERENCE=$nref OUTPUT=$dict2");
  system("java -jar $readGroups  INPUT=target2.newduped.bam OUTPUT=$out2 RGID=$lib RGLB=beads RGPL=illumina RGPU=lane2 RGSM=$lib");
  system("samtools index $out2");
system("rm target.duped.bam target.sorted.bam target.genome target.sorted.bam.bai target_pair.sam target_solo.sam outPairedSam1 outSoloSam1 target_pair.bam target_solo.bam target.metric target.bam target.newduped.bam target.newduped.bam.bai");

## run snpcaller again
    my $contigDir = $dir;
    my $nvcf = $contigDir .$lib . '.newFilter2.vcf';
    system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T HaplotypeCaller -R $newFasta -I $out2 -dt NONE -o $contigDir$lib.nhc.vcf");
    system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T ReadBackedPhasing -R $newFasta -I $out2 --variant $contigDir$lib.nhc.vcf --min_base_quality_score 10 -dt NONE -o $contigDir$lib.nraw.vcf");
    system("java -Djava.io.tmpdir=$tempDir -Xmx8g -jar $gatk -T VariantAnnotator -A DepthPerAlleleBySample -A AlleleBalance -A FisherStrand -A HaplotypeScore -A HardyWeinberg -R $newFasta -I $out2 --variant $contigDir$lib.nraw.vcf -o $contigDir$lib.nannotated.vcf");
    system("java -Djava.io.tmpdir=$tempDir -jar $gatk -T VariantFiltration -R $newFasta -o $nvcf --variant $contigDir$lib.nannotated.vcf $exp $filterQual $name \"qual\" $exp $filterDp20 $name \"dp20\" $exp $filterDp10 $name \"dp10\" $exp $filterHethigh $name \"hethigh\" $exp $filterHetlow $name \"hetlow\" $exp $filterAltpos $name \"altpos\" $exp $filterAltdef $name \"altdef\" $exp $filterHomhigh $name \"homhigh\" $exp $filterHomlow $name \"homlow\"");

### filter new vcf
    my $vcfv4 = $dir. $lib . '.short4.vcf';
    my $vcfv5 = $dir. $lib . '.short5.vcf';
    my $vcfv6 = $dir. $lib . '.short6.vcf';
##this ignores indels
    system("grep -v \'\#\' $nvcf \| awk \'\$8 ~ /ABH/\' > $vcfv4");

    open(VCF2, "<$vcfv4");
    open(NEWVCF2, ">>$vcfv5");
    foreach my $vcline2 (<VCF2>)    {
    my @tmp = split(/\t/,$vcline2);
    my $prea0; my $prea1; my $join; my $newletter; my $newline;
    my $a0; my $a1;
## for coverage less than 10, turn to Ns unless qual, then leave ref, or homozygous 1.0 then leave alt
if ($tmp[6] =~ /dp10/){
    if ($tmp[6] !~ /altdef|qual/){
$newletter = 'N';
$newline = $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$newletter."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\t".$tmp[8]."\t".$tmp[9]."\n";
print NEWVCF2 $newline;
    } #sub $tmp[4] to N
    if ($tmp[6] =~ m/altdef/) {
if ($tmp[6] !~ m/qual/) {
    print NEWVCF2 $vcline2;
}
    }
}
### for coverage less than 20, turn to IUPAC, using consensus later, unless qual, then leave ref, or homozygous 1.0, then leave alt
if ($tmp[6] =~ m/dp20/) {
    if ($tmp[6] !~ m/altdef|qual/){
$prea0 = $tmp[3];
$prea1 = $tmp[4];
$join = "$prea0$prea1";
unless (length($join) > 2) {
    for ($join ){
if (/TC|CT/) {
    $newletter = 'Y';
}  #sub $tmp[4] to Y
if (/AG|GA/){
    $newletter = 'R';
}  #sub $tmp[4] to R
if (/CA|AC/){
    $newletter = 'M';
}  #sub $tmp[4] to M
if (/TG|GT/){
    $newletter =  'K';
}  #sub $tmp[4] to K
if (/TA|AT/){
    $newletter = 'W';
} #sub $tmp[4] to W
if (/GC|CG/){
    $newletter = 'S';
} #sub $tmp[4] to S
$newline = $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$newletter."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\t".$tmp[8]."\t".$tmp[9]."\n";
print NEWVCF2 $newline;
    }
}
    }
    if ($tmp[6] =~  m/altdef/) {
print NEWVCF2 $vcline2;
    } #print
              
    if ($tmp[6] =~ m/qual/){
               
if ($tmp[6] =~ m/hetlow|altpos|homhigh|homlow|altdef/) {
    print NEWVCF2 $vcline2;
} #print
    }
}
 ## for qual, use 50 percent rule
               
if ($tmp[6] =~ m/qual/) {
                
    if ($tmp[6] !~ m/dp20|dp10/) {
if ($tmp[6] =~ m/hetlow|altpos|altdef|homhigh|homlow/) {
    print NEWVCF2 $vcline2;
} #print        
    }              
}
### no dp and no qual

if ($tmp[6] !~ m/qual|dp20|dp10/) {

    if ($tmp[6] =~ m/hetlow|homhigh|altdef/) {
print NEWVCF2 $vcline2;
    } #print
     
    if ($tmp[6] =~ m/altpos|PASS|homlow/) {
if ($tmp[7] !~ m/PhasingInconsistent/){
   $prea0 = $tmp[3];
    $prea1 = $tmp[4];
    $join = "$prea0$prea1";

    unless (length($join) > 2) { 
for ($join ){
    if (/TC|CT/) {
$newletter = 'Y';
    }  #sub $tmp[4] to Y
    if (/AG|GA/){
$newletter = 'R';
    }  #sub $tmp[4] to R
    if (/CA|AC/){
$newletter = 'M';
    }  #sub $tmp[4] to M
    if (/TG|GT/){
$newletter =  'K';
    }  #sub $tmp[4] to K
    if (/TA|AT/){
$newletter = 'W';
    } #sub $tmp[4] to W
    if (/GC|CG/){
$newletter = 'S';
    } #sub $tmp[4] to S
    $newline = $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$newletter."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\t".$tmp[8]."\t".$tmp[9]."\n";
   print NEWVCF2 $newline;
}
    }

}
                   
    }
} 


### deal with phasing inconsistent sites -- convert to N's, if want to count - count in original vcf with bed file

unless ($tmp[6] =~ m/qual|dp10|dp20/) { 
    for ($tmp[6] =~ m/PASS|homlow|altpos/) {
if ($tmp[7] =~ m/PhasingInconsistent/) {
   $newletter = 'N';
    $newline = $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\t".$tmp[3]."\t".$newletter."\t".$tmp[5]."\t".$tmp[6]."\t".$tmp[7]."\t".$tmp[8]."\t".$tmp[9]."\n";
    print NEWVCF2 $newline;
}   #sub $tmp[4] to D
    }
}
    }
    close VCF2;
    close NEWVCF2;

#####
### now use vcf-consensus to print out new file
 my $new2Fasta = $dir . $lib . '.newFasta2.fa';
   system("cat $columnName $vcfv5 | grep -v $getnoemptylines > $vcfv6");
   my $zipped2 = $vcfv6 . '.gz';
system("bgzip $vcfv6");
system("tabix -p vcf $zipped2");
system("cat $nref \| vcf-consensus $zipped2 > $new2Fasta");
system("samtools faidx $new2Fasta");
print "now finished with vcf - going to coverage\n";


## now make under 1x coverage bedfile and mask fasta for under 5x coverage, then cut out under 1x coverage from masked fasta                  

    my $tempcov0a = $dir . $lib . '.coverage0a';
    my $tempcov0b = $dir . $lib . '.coverage0b';
    my $tempcov5a = $dir . $lib . '.coverage5a';
    my $tempcov5b = $dir . $lib . '.coverage5b';
    my $over4fasta = $dir . $lib . '.over4.fasta';
    my $masked = $dir . $lib . '.covMask.fasta';
    my $seqlist = $dir . $lib . '.seqlist';
#make 0 cov bed file
    system("bedtools genomecov -ibam $out2 -d -split | awk \'\$3 > 0\' \| awk '\{print \$1\"\\t\"\(\$2-1\)\"\\t\"\$2\}' > $tempcov0a");
    system("bedtools merge -i $tempcov0a > $tempcov0b");
    print "made 0 coverage bed file\n";
#make less than 5 cov bedfile
    system("bedtools genomecov -ibam $out2 -d -split | awk \'\$3 < 5\' \| awk '\{print \$1\"\\t\"\(\$2-1\)\"\\t\"\$2\}' > $tempcov5a");
    system("bedtools merge -i $tempcov5a > $tempcov5b");
    print "made 5 coverage bed file\n";
#mask with less than 5x cov
system("bedtools maskfasta -fi $new2Fasta -bed $tempcov5b -fo $over4fasta");
  print "made masked 5x cov fasta\n";
#remove 0x cov
#select sequences from list of cov0b (these are sequences with any coverage at all)
    system("awk \'\{print \$1\}\' $tempcov0b | sort | uniq > $seqlist");
    system("$selectseqs -f $seqlist $new2Fasta > $masked");
    print "selected seqs with coverage to new file\n";
#system("rm $tempcov0a $Tempcov5a");
}
