#!/bin/bash
#SBATCH -J STAR
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH -o ~/Logs/%j.out
#SBATCH -e ~/Logs/%j.err
#SBATCH -n 12
##SBATCH --exclusive

echo "Mission starts at `date`" &&\

inputdir=~/projects/34_sh_rnaseq/01.RNA-seq
outdir=~/projects/34_sh_rnaseq/02.star2pass.out
genomeindex=~/software/referrence_files/index_star_2.7.8a_GRCh38_gencode.v22
annotation=~/software/referrence_files/gencode.v22.annotation.gtf

thread=12
readFileComm=zcat
R1=1.fq.gz
R2=2.fq.gz

cat ${1} |while read sample
do
STAR --genomeDir ${genomeindex} --readFilesCommand ${readFileComm} \
--readFilesIn ${inputdir}/${sample}_${R1} ${inputdir}/${sample}_${R2} \
--runThreadN ${thread} --alignIntronMax 500000 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 --alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory --limitBAMsortRAM 700000000000 \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 \
--outSAMstrandField intronMotif --outSAMattributes NH HI NM MD AS XS \
--sjdbGTFfile ${annotation} \
--limitSjdbInsertNsj 2000000 \
--outSAMmultNmax 1 \
--outSAMunmapped None \
--outFileNamePrefix ${outdir}/${sample}. \
--outBAMsortingThreadN ${thread} \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic 
done

echo "Mission ends at `date`" &&
