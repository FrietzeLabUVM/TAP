#!/bin/bash

#pools the input comma delimited list of bam files into the output file
#arg 1 is comma separated list of bam files
#arg 2 is name of output pooled bam file
if [ -z $INPUT ]; then
INPUT=$1
fi
if [ -z $OUTPUT ]; then
OUTPUT=$2
fi
bams=$INPUT
pooled=$OUTPUT
echo pooling $bams to $pooled
bams=${bams//","/" "}
topool=( $bams )
cmd="no command set"
if [ -f $pooled ]; then
        #echo pooled file already present, please verify $pooled is complete.
        cmd="echo $pooled already exists, check for completeness"
elif [ ${#topool[@]} -eq 1 ]; then
        echo pooling not necessary, just link for $bams to $pooled
        bams=$(readlink -f $bams)
        cmd="ln -s $(basename $bams) $pooled"
else
        echo gonna pool $bams into $pooled
        cmd="samtools merge $pooled $bams"
fi
echo CMD is "$cmd"
$cmd
for tp in ${topool[@]}; do
        echo $tp
done
echo pooling finished into $pooled
if [ -f "$pooled".bai ]; then
        echo indexing not necessary, already done
else
        echo indexing $pooled
        samtools index $pooled
fi
echo done

