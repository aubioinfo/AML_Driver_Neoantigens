#!/bin/bash
#SBATCH -J pVACfuse
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH -o ~/Logs/%j.out
#SBATCH -e ~/Logs/%j.err
#SBATCH -n 1
##SBATCH --exclusive

cache_VEP=~/.vep
path_FASTA=~/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
iedb=~/jinpeng/04.neoantigen/00_reAnalyze/data/iedb/iedb.fasta
maf2vcf=~/vcf2maf-1.6.18/maf2vcf.pl
path_VEP=~/ensembl-vep
path_VEP_Plugins=~/.vep/VEP_plugins
vep_run=${path_VEP}/vep
path_MHC=~/jinpeng/04.neoantigen/
assembly_version="GRCh37"
method="NetMHCpan"
epitope_len="8,9,10,11"
main_dir=~/jinpeng/04.neoantigen/driver_landscape
out_dir=${main_dir}/V5_driver_fusion10hla..pVACfuse.out
path_HLA=${main_dir}/scripts/10hla.txt

sample=sample
~/.conda/envs/pipeline/bin/pvacfuse run \
${main_dir}/V4_driver_fusion/${sample}/ \
${sample} \
$(grep ${sample} ${path_HLA} | awk '{print $2}') \
${method} \
${out_dir}/${sample} \
-e1 ${epitope_len} \
-t 4 \
-k \
--iedb-install-directory ${path_MHC}