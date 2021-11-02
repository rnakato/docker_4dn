#!/bin/bash
coolfile=$1
outdir=$2
binsize=$3
gt=$4

tmpfile=$(mktemp "/tmp/${0##*/}.tmp.XXXXXX")

pwd=$(cd $(dirname $0) && pwd)
chrlist=$($pwd/getchr_from_genometable.sh $gt)

dir=$outdir/Matrix/$binsize
mkdir -p $dir
for chr in $chrlist
do
    if test $chr = "chrY" -o $chr = "chrM" -o $chr = "chrMT" ;then continue; fi

    echo $chr
    cooler dump -b --na-rep NA --join -r $chr $coolfile::resolutions/$binsize | cut -f2,5,8 | grep -v NA > $tmpfile
    python3 $pwd/convert_JuicerDump_to_dense.py $tmpfile $dir/$chr.matrix.gz $gt $chr $binsize
    rm $tmpfile
done
