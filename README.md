# docker_4dn
Repository of Docker image for 4dn pipeline. Docker image is at: https://hub.docker.com/repository/docker/rnakato/4dn/general

## Run

For Docker:

    # pull docker image
    docker pull rnakato/4dn

    # container login
    docker run [--gpus all] --rm -it rnakato/4dn /bin/bash
    # jupyter notebook
    docker run [--gpus all] --rm -p 8888:8888 -v (your directory):/opt/work rnakato/4dn <command>

For Singularity:

    # build image
    singularity build rnakato_4dn.sif docker://rnakato/4dn 
    # jupyter notebook
    singularity exec [--nv] rnakato_4dn.sif <command>

## Usage

These scripts assume that the fastq files are stored in `fastq/$cell` (e.g., `fastq/Control_1`).
The outputs are stored in `Results4dn/$cell`.

### Mapping reads

    odir=Results4dn/$cell
    bwa_index=bwa-indexes/UCSC-hg38  # specify BWA index
    fq1s=""
    fq2s=""
    for id in $(ls fastq/$cell/*_R1.fastq.gz)
    do
        id=${id%%_R1.fastq.gz}
        fq1s="$fq1s ${id}_R1.fastq.gz"
        fq2s="$fq2s ${id}_R2.fastq.gz"
    done
    
    mkdir -p $odir/aligned/ $odir/log
    bwa mem -t 32 -SP5M $bwa_index <(zcat $fq1s) <(zcat $fq2s) | samtools view -Shb - > $odir/aligned/$cell.bam

### pairtools parse/sort/markdup

    odir=Results4dn/$cell
    gt=genome_table.txt   # specify genome_table
    BAM=$odir/aligned/$cell.bam
    OUTDIR=$odir/pairsam-parse-sort
    OUTPREFIX=$cell
    SORTED_PAIRS_PATH=${OUTDIR}/${OUTPREFIX}.sam.pairs.gz
    PAIRSAM=${OUTDIR}/${OUTPREFIX}.sam.pairs.gz
    MARKED_PAIRSAM=$OUTDIR/${OUTPREFIX}.marked.sam.pairs.gz
    sing="singularity exec rnakato_4dn.sif"

    mkdir -p $OUTDIR
    samtools view -h $BAM | {
        $sing pairtools parse -c $gt --add-columns mapq
    } | {
        $sing pairtools sort --nproc 32 \
              --tmpdir ${OUTDIR} \
              --output ${SORTED_PAIRS_PATH}
    }
    $sing pairtools dedup --mark-dups --output-dups - --output-unmapped - --output ${MARKED_PAIRSAM} ${PAIRSAM}
    $sing pairix ${MARKED_PAIRSAM}  # sanity check

### Filtering (select, split)

    odir=Results4dn/$cell
    OUTDIR=$odir/pairsam-parse-sort
    OUTPREFIX=$cell
    sing="singularity exec rnakato_4dn.sif"

    PAIRSAM=${OUTDIR}/${OUTPREFIX}.sam.pairs.gz
    UNMAPPED_SAMPAIRS=$OUTDIR/${OUTPREFIX}.unmapped.sam.pairs.gz
    DEDUP_PAIRS=$OUTDIR/${OUTPREFIX}.dedup.pairs.gz
    LOSSLESS_BAM=$OUTDIR/${OUTPREFIX}.lossless.bam
    TEMPFILE=$OUTDIR/temp.gz
    TEMPFILE1=$OUTDIR/temp1.gz

    ## Generate lossless bam
    $sing pairtools split --output-sam ${LOSSLESS_BAM} ${PAIRSAM}

    # Select UU, UR, RU reads
    $sing pairtools select '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
                 --output-rest ${UNMAPPED_SAMPAIRS} \
                 --output ${TEMPFILE} \
                ${PAIRSAM}
    $sing pairtools split --output-pairs ${TEMPFILE1} ${TEMPFILE}
    $sing pairtools select 'True' --chrom-subset $gt -o ${DEDUP_PAIRS} ${TEMPFILE1}
    $sing pairix ${DEDUP_PAIRS}
    rm ${TEMPFILE} ${TEMPFILE1}

### QC report generater

    odir=Results4dn/$cell
    OUTDIR=$odir/pairsam-parse-sort
    OUTPREFIX=$cell
    input_pairs=${OUTDIR}/${OUTPREFIX}.marked.sam.pairs.gz
    sample_name=$cell
    sing="singularity exec rnakato_4dn.sif"
    scriptdir=/usr/local/bin/pairsqc/
    
    $sing python3 $scriptdir/pairsqc.py -p $input_pairs -c $gt -tP -s $sample_name -O $odir/$sample_name -M $max_distance
    $sing Rscript $scriptdir/plot.r $enzymelen $odir/$sample_name\_report

### addfrag2pairs

    odir=Results4dn/$cell
    sing="singularity exec rnakato_4dn.sif"
    scriptdir=/usr/local/bin
    OUTDIR=$odir/pairsam-parse-sort
    OUTPREFIX=$cell
    input_pairs=${OUTDIR}/${OUTPREFIX}.marked.sam.pairs.gz
    gunzip -c $input_pairs \
        | $sing $scriptdir/pairix/util/fragment_4dnpairs.pl -a - ${OUTDIR}/${OUTPREFIX}.ff.pairs $restrictionsite
    $sing bgzip  -f ${OUTDIR}/${OUTPREFIX}.ff.pairs
    $sing pairix -f ${OUTDIR}/${OUTPREFIX}.ff.pairs.gz

## Cooler

    odir=Results4dn/$cell
    sing="singularity exec rnakato_4dn.sif"
    OUTDIR=$odir/coolfile
    OUTPREFIX=$cell
    pairs_file=$odir/pairsam-parse-sort/${OUTPREFIX}.ff.pairs.gz
    out_prefix=${OUTDIR}/${OUTPREFIX}
    max_split=2
    binsizes="5000,10000,25000,50000,100000,500000,1000000,2500000,5000000,10000000"
    binsize_min=5000
    pwd=`pwd`

    mkdir -p $OUTDIR
    $sing cooler cload pairix -p $ncore -s $max_split $gt:$binsize_min $pairs_file $out_prefix.cool " >& $odir/log/run-cooler.$binsize_min
    ### matrix balancing
    $sing cooler balance -p $ncore $out_prefix.cool
    $sing run-juicebox-pre.sh -i $pairs_file -c $gt -o $odir/$cell -r $binsize_min -u $binsizes
    ### Generate a multi-resolution cooler file by coarsening
    $sing run-cool2multirescool.sh -i $out_prefix.cool -p $ncore -o $out_prefix -u $binsizes

    ### Adds a normalization vector from a hic file to an mcool file
    $sing run-add-hicnormvector-to-mcool.sh $pwd/$odir/$cell.hic $pwd/$out_prefix.multires.cool $pwd/$OUTDIR

    ### dump matrix files
    for binsize in 5000 25000 50000 100000
    do
        cool=$out_prefix.$binsize.cool
        $sing cooler cload pairix -p $ncore -s $max_split $gt:$binsize $pairs_file $cool
        ### matrix balancing
        $sing cooler balance -p $ncore $cool
        $sing run-cooler-dump.sh $cool $odir $binsize $gt
    done

    cp $pwd/$out_prefix.multires.cool $out_prefix.hic2cool.cool
    $sing hic2cool extract-norms -e $odir/$cell.hic $out_prefix.hic2cool.cool

## Build image from Dockerfile

First clone and move to the repository

    git clone https://github.com/rnakato/docker_4dn.git
    cd docker_4dn

Then type:

    docker build -t <account>/4dn .

## Contact

Ryuichiro Nakato: rnakato AT iqb.u-tokyo.ac.jp
