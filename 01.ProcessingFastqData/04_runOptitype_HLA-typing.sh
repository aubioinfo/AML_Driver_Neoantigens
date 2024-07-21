### HLA-typing using optitype
# conda activate hlatyping
# module load miniconda3
# pip install future
# pip install matplotlib
# pip install pyomo
# pip install tables
left_file=R1.fq.gz
right_file=R2.fq.gz
threads=15
sample_file=tcga.list
input_dir=~/AML/beatAML/Fastq
output_dir1=~/jinpeng/neoantigen/00_Analyze/01.hla.typing/bam_fastq
output_dir2=~/jinpeng/neoantigen/00_Analyze/02.hla.typing/
refe_fasta=~/jinpeng/software/OptiType/data/hla_reference_rna.fasta
optitype_dir=~/.conda/envs/hlatyping/bin/

cat ${sample_file} |while read Sample
do
# Mapping the reads to the HLA reference
razers3 -i 95 -m 1 -dr 0 -tc ${threads} -o ${output_dir1}/${Sample}_R1.bam ${refe_fasta} ${input_dir}/${Sample}_${left_file}
razers3 -i 95 -m 1 -dr 0 -tc ${threads} -o ${output_dir1}/${Sample}_R2.bam ${refe_fasta} ${input_dir}/${Sample}_${right_file}
#Transform bam files to Fastq files
samtools bam2fq ${output_dir1}/${Sample}_R1.bam > ${output_dir1}/${Sample}_HLA_R1.fastq
samtools bam2fq ${output_dir1}/${Sample}_R2.bam > ${output_dir1}/${Sample}_HLA_R2.fastq
# Run optitype
python ${optitype_dir}/OptiTypePipeline.py -i ${output_dir1}/${Sample}_HLA_R1.fastq ${output_dir1}/${Sample}_HLA_R2.fastq --rna -v -o ${output_dir2}/${Sample}
done

# Merging multiple-samples
ls `pwd` > all.dir
grep 'BA' all.dir > beat.list
grep 'TCGA' all.dir > tcga.list
cat beat.list tcga.list > list 
rm beat.list tcga.list

for Sample in `cat list`
do
	let i+=1
    printf "\rprocessing %s .." ${Sample}
    
	if [ ${i} -eq 1 ]; then
        cat ${Sample}/2021*/*.tsv > HLA.combined.tsv
    else
        cat ${Sample}/2021*/*.tsv | sed '1d' >> HLA.combined.tsv
    fi
done

awk -F "\t" '{print $2,$3,$4,$5,$6,$7}' HLA.combined.tsv > HLA.combined.list.txt
sed -i 's/A/HLA-A/g' HLA.combined.list.txt
sed -i 's/B/HLA-B/g' HLA.combined.list.txt
sed -i 's/C/HLA-C/g' HLA.combined.list.txt