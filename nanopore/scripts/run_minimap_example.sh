# Run minimap2 and convert to sorted BAM
minimap2 -ax map-ont -t 8 ./ref/GRCh38.primary_assembly.canonical.EBV.fa "${FASTQ_FILE}" | \
samtools sort -@ 8 -o "${OUT_BAM}"

samtools index "${OUT_BAM}"
