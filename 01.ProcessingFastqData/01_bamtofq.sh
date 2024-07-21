#!/bin/bash
NAME=$1
echo -n $NAME ...
samtools sort -n -@ 8 ~/upload/${NAME}.bam | samtools fastq -@ 8 -1 ~/AML/${NAME}_R1.fq.gz -2 /lustre/home/acct-medkkw/medzw/AML/${NAME}_R2.fq.gz -0 /dev/null -s /dev/null -
echo Done.