#!/bin/bash
#SBATCH -J htseq
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH -o ~/Logs/%j.out
#SBATCH -e ~/Logs/%j.err
#SBATCH -n 1
##SBATCH --exclusive

echo "Mission starts at `date`" &&\

# conda activate base
annotation=~/referrence_sorf.index/gencode.v22.annotation.gtf
inputdir=~/projects/34_sh_rnaseq/02.star2pass.out
outdir=~/projects/34_sh_rnaseq/03.htseq.out

sample=$1

htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty ${inputdir}/${sample}.Aligned.sortedByCoord.out.bam ${annotation} > ${outdir}/${sample}.count

echo "Mission ends at `date`" 