#!/bin/bash
#SBATCH -J pvacTo
#SBATCH -p small
#SBATCH --mail-type=end
#SBATCH --mail-user=
#SBATCH -o ~/Logs/%j.out
#SBATCH -e ~/Logs/%j.err
#SBATCH -n 4
##SBATCH --exclusive

cache_VEP=~/.vep
path_FASTA=~/.vep/homo_sapiens/101_GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
iedb=~/04.neoantigen/00_reAnalyze/data/iedb/iedb.fasta
maf2vcf=~/vcf2maf-1.6.18/maf2vcf.pl
path_VEP=~/ensembl-vep
path_VEP_Plugins=~/.vep/VEP_plugins
vep_run=${path_VEP}/vep
path_MHC=~/04.neoantigen/
assembly_version="GRCh37"
method="NetMHCpan"
epitope_len="8,9,10,11"
main_dir=~/04.neoantigen/driver_landscape
maf_file=${main_dir}/scripts/hg19_example.vep.maf
out_dir=~/04.neoantigen/driver_landscape/V3.aml.Driver
path_HLA=${main_dir}/scripts/aml_patients.hla 
VEP_annotated_vcf=${out_dir}/VEP_annotated_vcf
pvacseq_out=${out_dir}/pvacseq_Results_netmhcPan

# Check if the out_dir exists
if [ ! -d "${out_dir}" ]; then
    echo  -e "\033[31m The output directory does not exist, creating it ! \033[0m" 
    mkdir -p ${out_dir}
fi
# Check if the tumor-normal directory exists
tumor_normal=${out_dir}/tumor_normal_pair
if [ ! -d "${tumor_normal}" ]; then
    echo  -e "\033[31m The tumor-normal-pair directory not exist under ${out_dir}, creating it ! \033[0m" 
    mkdir -p ${tumor_normal}
fi

# transform
echo "[------Transforming started--------]"
perl ${maf2vcf} --input-maf ${maf_file} --output-dir ${tumor_normal} --ref-fasta ${path_FASTA} --per-tn-vcfs
echo "[------Transforming completed------]"

# deleting the normal-sample column
echo "[------Deleting the column of normal-sample------]"
single_vcf=${out_dir}/single_vcfs
if [ ! -d "${sgvcf}" ]; then
    echo  -e "The single_vcfs directory does not exist under ${out_dir}, creating it !" 
    mkdir -p ${single_vcf}
fi

for vcf in $(ls ${tumor_normal} | grep "vs")
do
    Sample=$(echo ${vcf} | sed -E 's/^(.*)_vs.*/\1/')
    printf "\rprocessing %s .." ${Sample}
    cat ${tumor_normal}/$vcf | awk 'BEGIN{FS="\t";OFS="\t"}{if($1 ~ /^##]/){print $0}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}}' > ${single_vcf}/${Sample}".vcf"
done

echo 
echo "[------Transforming step finished successfully !------]"

# annotating VCF files using VEP
VEP_annotated_vcf=${out_dir}/VEP_annotated_vcf
if [ ! -d "${VEP_annotated_vcf}" ]; then
    echo  -e "The VEP_annotated_vcf directory does not exist under ${out_dir}, creating it !" 
    mkdir -p ${VEP_annotated_vcf}
fi
# obtaining all input vcf filenames
filenames=$(ls "${single_vcf}" | grep -E '.*\.vcf$')

for filename in ${filenames}
do
sampleID=$(echo ${filename} | sed -E 's/(.*)\.vcf$/\1/')
printf "\rprocessing %s ..\n" ${sampleID}
${vep_run} --input_file ${single_vcf}/${filename} \
--format vcf --output_file stdout \
--vcf --symbol --terms SO \
--plugin Frameshift --plugin Wildtype \
--dir_plugins ${path_VEP_Plugins} \
--assembly ${assembly_version} \
--fasta ${path_FASTA} \
--dir_cache ${cache_VEP} \
--offline --cache_version 101 \
--pick --force_overwrite \
> ${VEP_annotated_vcf}/${sampleID}"_VEP_annotated_filterd.vcf"
done
# runing pvacseq to predict the neoantgens
pvacseq_out=${out_dir}/pvacseq_Results_netmhcPan
if [ ! -d "${pvacseq_out}" ]; then
    echo  -e "\033[33m The pvacseq_out directory does not exist under ${out_dir}, creating it. \033[0m" 
    mkdir -p ${pvacseq_out}
fi 

filenames=$(ls "${single_vcf}" | grep -E '.*\.vcf$')
for filename in ${filenames}
do
sampleID=$(echo ${filename} | sed -E 's/(.*)\.vcf$/\1/')
printf "\rprocessing %s ..\n" ${sampleID}
~/.conda/envs/pipeline/bin/pvacseq run \
${VEP_annotated_vcf}/${sampleID}"_VEP_annotated_filterd.vcf" \
${sampleID} \
$(grep $sampleID ${path_HLA} | awk '{print $2}') \
NetMHCpan \
${pvacseq_out}/${sampleID} \
--normal-cov 0 \
--tdna-cov 0 \
--trna-cov 0 \
--normal-vaf 0 \
--tdna-vaf 0 \
--trna-vaf 0 \
--expn-val 0 \
-e1 $epitope_len \
-a sample_name \
-b 500 \
-k \
-t 20 \
--iedb-install-directory ${path_MHC}
done