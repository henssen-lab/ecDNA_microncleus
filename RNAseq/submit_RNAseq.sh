sample=(RNA_01_DM_HU_1_S1 RNA_02_DM_HU_2_S2 RNA_03_DM_HU_3_S3  RNA_04_DM_UT_1_S4 RNA_05_DM_UT_2_S5 RNA_06_DM_UT_3_S6 )


trim_galore --nextseq 20 --phred33 --illumina -j 8 -o /path/to/output_dir /path/to/fastq_dir

STAR --genomeDir /path/to/STAR --readFilesIn ${sample[${SLURM_ARRAY_TASK_ID}]}_trimmed.fq.gz --runThreadN 64 --genomeLoad NoSharedMemory --readFilesCommand zcat --twopassMode Basic --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sample[${SLURM_ARRAY_TASK_ID}]}

htseq-count -r pos -s no ${sample[${SLURM_ARRAY_TASK_ID}]}Aligned.sortedByCoord.out.bam /data/cephfs-1/work/projects/ecdna-chrom/reference/GRCh38_canonical_EBV/gencode.v47.primary_assembly.annotation.gtf > ${sample[${SLURM_ARRAY_TASK_ID}]}.htseq.counts


