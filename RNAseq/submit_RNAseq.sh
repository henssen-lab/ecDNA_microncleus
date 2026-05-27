#!/bin/bash

#SBATCH --job-name=RNA_hepato_PDX
#SBATCH --ntasks=64
#SBATCH --nodes=1
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=robin.xu@charite.de
#SBATCH --mail-type=ALL
#SBATCH --output=/data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/scripts/log/R-%x.%j.out
#SBATCH --error=/data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/scripts/log/R-%x.%j.err



eval "$($(which conda) shell.bash hook)"
conda activate chipseq
set -x



sample=(RNA_01_DM_HU_1_S1 RNA_02_DM_HU_2_S2 RNA_03_DM_HU_3_S3  RNA_04_DM_UT_1_S4 RNA_05_DM_UT_2_S5 RNA_06_DM_UT_3_S6 )



#cd /data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/raw_data/A4713/250109_LH00253_0191_A225K25LT1/

#cat A4713_${sample[${SLURM_ARRAY_TASK_ID}]}_L001_R1_001.fastq.gz A4713_${sample[${SLURM_ARRAY_TASK_ID}]}_L002_R1_001.fastq.gz > ${sample[${SLURM_ARRAY_TASK_ID}]}.fastq.gz



#mkdir /data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/raw_data/A4713/out/data/${sample[${SLURM_ARRAY_TASK_ID}]}
cd /data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/raw_data/A4713/out/data/${sample[${SLURM_ARRAY_TASK_ID}]}


#trim_galore --nextseq 20 --phred33 --illumina -j 8 -o /data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/raw_data/A4713/out/data/${sample[${SLURM_ARRAY_TASK_ID}]} /data/cephfs-1/work/projects/ecdna-chrom/RNAseq/MN_RNA/raw_data/A4713/250109_LH00253_0191_A225K25LT1/${sample[${SLURM_ARRAY_TASK_ID}]}.fastq.gz

#STAR --genomeDir /data/cephfs-1/work/projects/ecdna-chrom/reference/GRCh38_canonical_EBV/STAR/ --readFilesIn ${sample[${SLURM_ARRAY_TASK_ID}]}_trimmed.fq.gz --runThreadN 64 --genomeLoad NoSharedMemory --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample[${SLURM_ARRAY_TASK_ID}]}

htseq-count -r pos -s no ${sample[${SLURM_ARRAY_TASK_ID}]}Aligned.sortedByCoord.out.bam /data/cephfs-1/work/projects/ecdna-chrom/reference/GRCh38_canonical_EBV/gencode.v47.primary_assembly.annotation.gtf > ${sample[${SLURM_ARRAY_TASK_ID}]}.htseq.counts


