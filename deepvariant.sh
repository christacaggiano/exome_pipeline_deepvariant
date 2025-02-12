
#!/bin/bash
#$ -cwd
#$ -r y
#$ -j y
#$ -l h_data=80G
#$ -l h_rt=24:00:00
#$ -l highp
#$ -t 2-11

######

########################### environment #####################################


source ~/.bashrc
conda activate /u/home/c/christac/project-nzaitlen/miniconda/envs/exome

ID_NUM=$((SGE_TASK_ID+40))

# Input parameters
FASTQ_R1="/u/project/zaitlenlab/christac/exome_sequencing/trimmed/"$ID_NUM"_R1_val_1.fq.gz" 
FASTQ_R2="/u/project/zaitlenlab/christac/exome_sequencing/trimmed/"$ID_NUM"_R2_val_2.fq.gz" 
REFERENCE="/u/project/zaitlenlab/christac/exome_sequencing/t2t/chm13v2.0.fa"
OUTPUT_DIR="/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2" 
SAMPLE_NAME="sample"$ID_NUM 
THREADS=48 
DEEPVARIANT_SIF="/u/project/zaitlenlab/christac/exome_sequencing/sif/deepvariant.sif" 
TMPDIR="/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/temp/"$ID_NUM 
GLNEXUS_SIF="/u/project/zaitlenlab/christac/exome_sequencing/sif/glnexus.sif"

########## Step 1: Crleate output directory
mkdir -p ${OUTPUT_DIR}

######### Step 2: Index the reference genome (if not already indexed)
if [ ! -f ${REFERENCE}.bwt ]; then
    echo "Indexing reference genome..."
    bwa index ${REFERENCE}
fi

########## Step 3: Align reads with BWA
echo "Aligning reads with BWA..."

echo ${FASTQ_R1} 
echo ${FASTQ_R2}

bwa mem -t ${THREADS} ${REFERENCE} ${FASTQ_R1} ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE_NAME}.sam

########## Step 4: Convert SAM to BAM and sort
echo "Converting SAM to BAM and sorting..."
samtools view -Sb ${OUTPUT_DIR}/${SAMPLE_NAME}.sam > ${OUTPUT_DIR}/${SAMPLE_NAME}.bam
samtools sort -o ${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam ${OUTPUT_DIR}/${SAMPLE_NAME}.bam

rm ${OUTPUT_DIR}/${SAMPLE_NAME}.sam 
rm ${OUTPUT_DIR}/${SAMPLE_NAME}.bam

########### Step 5: Mark duplicates with Picard
echo "Marking duplicates with Picard..."
picard MarkDuplicates \
    I=${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam \
    O=${OUTPUT_DIR}/${SAMPLE_NAME}_marked_duplicates.bam \
    M=${OUTPUT_DIR}/${SAMPLE_NAME}_marked_dup_metrics.txt

# ########### Step 6: Index the BAM file
echo "Indexing the BAM file..."
samtools index ${OUTPUT_DIR}/${SAMPLE_NAME}_marked_duplicates.bam

########### Step 7: Run DeepVariant using Docker
echo "Running DeepVariant..."

rm -rf ${TMPDIR}
mkdir -p ${TMPDIR}

singularity exec \
    --bind  /u/project/zaitlenlab/christac/exome_sequencing,:/mnt,${TMPDIR}:/tmp \
    ${DEEPVARIANT_SIF} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref /u/project/zaitlenlab/christac/exome_sequencing/t2t/chm13v2.0.fa   \
    --reads "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/"$SAMPLE_NAME"_marked_duplicates.bam" \
    --output_vcf "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/"$SAMPLE_NAME".vcf.gz" \
    --output_gvcf "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/"$SAMPLE_NAME".g.vcf" \
    --intermediate_results_dir "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/temp/"$ID_NUM \
    --num_shards=${THREADS}

############  Step 8: Index the VCF file
echo "Indexing the VCF file..."
tabix -f vcf ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf.gz

############  Step 9: Quality control statistics 
bcftools stats ${OUTPUT_DIR}/${SAMPLE_NAME}.vcf.gz > ${OUTPUT_DIR}/${SAMPLE_NAME}_stats.txt

rm -rf ${TMPDIR}
mkdir -p ${TMPDIR}

singularity exec \
    --bind  /u/project/zaitlenlab/christac/exome_sequencing,:/mnt,${TMPDIR}:/tmp \
    --env TMPDIR="/tmp" \
    ${DEEPVARIANT_SIF} \
    /opt/deepvariant/bin/vcf_stats_report \
    --input_vcf "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/"$SAMPLE_NAME".vcf.gz" \
    --outfile_base "/u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2/"$SAMPLE_NAME 
    

############ Step 10: Rename and merge

old_sample="default"
new_sample=${SAMPLE_NAME}

sed "s/${old_sample}/${new_sample}/" deepvariant_v2/${SAMPLE_NAME}.g.vcf > deepvariant_v2/${SAMPLE_NAME}_name.g.vcf

singularity shell --bind /u/project/zaitlenlab/christac/exome_sequencing/deepvariant_v2:/data sif/glnexus.sif
glnexus_cli --config DeepVariantWES /data/*.g.vcf > combined.bcf














