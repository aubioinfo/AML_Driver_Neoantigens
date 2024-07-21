#!/bin/bash
#SBATCH -J maf2maf
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=1442518429@qq.com
#SBATCH -o ~/Logs/%j.out
#SBATCH -e ~/Logs/%j.err
#SBATCH -n 4
##SBATCH --exclusive

echo "Mission starts at `date`" &&\

perl /lustre/home/acct-medkkw/medzw/vcf2maf-1.6.18/maf2maf.pl --input-maf AML_driverGenes.maf --output-maf AML_driverGenes.VEP.maf --ref-fasta ~/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa --vep-path=$HOME/ensembl-vep 

echo "Mission completes at `date`"