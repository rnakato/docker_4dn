#!/bin/bash

ex(){ echo $1; eval $1; }

fq1s=""
fq2s=""

for cell in $(ls fastq/* -d); do
cell=$(basename $cell)
for id in $(ls fastq/$cell/*_1.fastq.gz)
do
    id=${id%%_1.fastq.gz}
#    echo $id
    fq1s="$fq1s ${id}_1.fastq.gz"
    fq2s="$fq2s ${id}_2.fastq.gz"
done

bwa_index=/work/Database/bwa-indexes/UCSC-hg38
ncore=24
odir=Results4dn/$cell

ex "mkdir -p $odir/aligned/"
ex "bwa mem -t $ncore -SP5M $bwa_index <(zcat $fq1s) <(zcat $fq2s) | samtools view -Shb - | samtools sort > $odir/aligned/$cell.sorted.bam"
ex "samtools index $odir/aligned/$cell.sorted.bam"

singularity="singularity exec /work/SingularityImage/rnakato_4dn.img"

$singularity run-bam2pairs.sh $odir/aligned/$cell.sorted.bam $odir/aligned/$cell

pairfile=$odir/aligned/$cell.bsorted.pairs.gz
gt=/work/Database/UCSC/hg38/genome_table
restrictionsite=/work/Database/HiC-restriction_sites/MboI_resfrag_hg38.bed
resolution=100000

$singularity run-cooler.sh $pairfile $gt $resolution $ncore $odir/$cell 2

$singularity run-cool2multirescool.sh -i $odir/$cell.cool -p $ncore -o $odir/$cell
$singularity run-pairsqc-single.sh $pairfile $gt $cell 4 $odir/
$singularity run-addfrag2pairs.sh $pairfile $restrictionsite $odir/aligned/$cell
$singularity run-juicebox-pre.sh -i $pairfile -c $gt -o $odir/$cell
$singularity run-add-hicnormvector-to-mcool.sh `pwd`/$odir/$cell.hic $odir/$cell.multires.cool $odir
cp $odir/$cell.multires.cool $odir/$cell.mcool
$singularity hic2cool extract-norms -e $odir/$cell.hic $odir/$cell.mcool

$singularity cooler show --out test.png --dpi 200 Results4dn/MCF-7_E2minus_rep1/MCF-7_E2minus_rep1.cool chr3:0-80,000,000
done
