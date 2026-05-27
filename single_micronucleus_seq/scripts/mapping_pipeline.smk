#Snakemake workflow for trimming and mapping single cell data. Aligns to hg38.

SAMPLES = []
for root, dirs, files in os.walk("./sMNseq/pilot/data/", topdown=False):
    for name in files:
        if os.path.join(root, name).endswith('_L001_R1_001.fastq.gz'):
            SAMPLES.append(os.path.join(root, name).removesuffix('_L001_R1_001.fastq.gz').removeprefix('./sMNseq/pilot/data/'))

print(SAMPLES)
rule all:
    input:
        expand('./sMNseq/pilot/samples/{sample}/{sample}_L001_R1_001_val_1.fq.gz', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/{sample}_L001_R2_001_val_2.fq.gz', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam.bai', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam.bai', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_metrics.txt', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam.bai', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam.bw', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_CPM.bam.bw', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_CPM.bam.bedgraph', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}.bed', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_circ_counts.bed', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_lin_counts.bed', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_lin_circ_counts.txt', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_25kb_bin.bed', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_q20_dedup.bam', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_25kb_bin_q20.bed', sample=SAMPLES),
        expand('./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_unique_read.count', sample=SAMPLES)


rule trim:
    threads: 8
    resources:
        mem='64G',
        time='4:00:00',
    input:
        read1='./sMNseq/pilot/data/{sample}_L001_R1_001.fastq.gz',
        read2='./sMNseq/pilot/data/{sample}_L001_R2_001.fastq.gz',
    params:
        out_dir='./sMNseq/pilot/samples/{sample}/'
    output:
        read1='./sMNseq/pilot/samples/{sample}/{sample}_L001_R1_001_val_1.fq.gz',
        read2='./sMNseq/pilot/samples/{sample}/{sample}_L001_R2_001_val_2.fq.gz'
    shell:
        'trim_galore --nextseq 20 --phred33 --illumina --paired -j 8 -o {params.out_dir} {input.read1} {input.read2}'

rule bwa:
    input:
        read1='./sMNseq/pilot/samples/{sample}/{sample}_L001_R1_001_val_1.fq.gz',
        read2='./sMNseq/pilot/samples/{sample}/{sample}_L001_R2_001_val_2.fq.gz',
        genome= '/data/cephfs-1/work/projects/ecdna-chrom/reference/GRCh38_canonical_EBV/GRCh38.primary_assembly.canonical.EBV.fa'
    output:
        sam=temp('./sMNseq/pilot/samples/{sample}/{sample}.sam')
    resources:
        mem="64G",
        time="4:00:00",
    threads: 16
    log: "./sMNseq/pilot/log/{sample}_bwamem.log"
    shell:
        'bwa mem -q -t {threads} {input.genome} {input.read1} {input.read2} > {output.sam} 2>{log}'

rule coordinate_sort:
    input:
        sam='./sMNseq/pilot/samples/{sample}/{sample}.sam'
    output:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam'
    resources:
        mem="32G",
        time="4:00:00",
    threads: 8
    log: "./sMNseq/pilot/log/{sample}_sort.log"
    shell:
        'samtools sort -@ {threads} -o  {output.sorted_bam} {input.sam} 2>{log}'

rule index_sorted_bam:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam'
    output:
        sorted_bam_idx='./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam.bai'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'samtools index {input.sorted_bam} 2>{log}'

rule add_rg:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}.bam'
    output:
        sorted_bam_rg='./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam'
    params:
        name="{sample}"
    log: "./sMNseq/pilot/log/{sample}_addrg_sorted.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'samtools addreplacerg -r "@RG\tID:RG1\tSM:{params.name}\tPL:Illumina" -o {output.sorted_bam_rg} {input.sorted_bam} 2>{log}'


rule index_sorted_bam_rg:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam'
    output:
        sorted_bam_idx='./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam.bai'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted_rg.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'samtools index {input.sorted_bam} 2>{log}'

rule mark_dup:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam',
        sorted_bam_idx='./sMNseq/pilot/samples/{sample}/sorted_{sample}_rg.bam.bai'
    output:
        deduped_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam',
        dup_metrics='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_metrics.txt'
    log: "./sMNseq/pilot/log/{sample}_dedup.log"
    resources:
        mem="64G",
        time="4:00:00",
    threads: 8
    shell:
        'picard MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT={input.sorted_bam} OUTPUT={output.deduped_bam} ASSUME_SORTED=true METRICS_FILE={output.dup_metrics} REMOVE_DUPLICATES=true 2>{log}'

rule index_sorted_bam_dedup:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam'
    output:
        sorted_bam_idx='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam.bai'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted_dedup.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'samtools index {input.sorted_bam} 2>{log}'


rule coverage:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam'
    output:
        sorted_bam_covearge='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam.bw'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted_dedup.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'bamCoverage -b {input.sorted_bam} -o {output.sorted_bam_covearge}  2>{log}'

rule coverage_CPM:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam'
    output:
        sorted_bam_covearge='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_CPM.bam.bw'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted_dedup_CPM.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'bamCoverage -b {input.sorted_bam} -bs 50000 --normalizeUsing CPM -o {output.sorted_bam_covearge}  2>{log}'

rule coverage_bedgraph:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam'
    output:
        sorted_bam_covearge='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_CPM.bam.bedgraph'
    log: "./sMNseq/pilot/log/{sample}_idx_sorted_dedup_CPM_bedgraph.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        'bamCoverage -b {input.sorted_bam} -bs 1000000 --normalizeUsing CPM --minMappingQuality 20 --centerReads -of bedgraph -o {output.sorted_bam_covearge}  2>{log}'


rule bin_25kb:
    input:
        bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam',
    output:
        bed='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_25kb_bin.bed'
    params:
        bed="./sMNseq/pilot/scripts/ref/GRCh38.primary_assembly.canonical.EBV.25kb_bins.bed"
    log: "./sMNseq/pilot/log/{sample}_count_lin.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        "bedtools coverage -a {params.bed} -b {input.bam} > {output.bed} 2>{log}"

rule bamtobam:
    input:
        sorted_bam='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup.bam'
    output:
        bam='./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_q20_dedup.bam'

    log: "./sMNseq/pilot/log/{sample}_bamtobam.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        "samtools view -b -q 20 -F 3844 {input.sorted_bam} > {output.bam} 2>{log}"

rule bin_25kb_bamtobam:
    input:
        bam='./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_q20_dedup.bam',
    output:
        bed='./sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_25kb_bin_q20.bed'
    params:
        bed="./sMNseq/pilot/scripts/ref/GRCh38.primary_assembly.canonical.EBV.25kb_bins.bed"
    log: "./sMNseq/pilot/log/{sample}_count_lin_bin_25kb_bamtobam.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        "bedtools coverage -a {params.bed} -b {input.bam} > {output.bed} 2>{log}"

rule count:
    input:
        bam='./sMNseq/pilot/samples/{sample}/bed/sorted_{sample}_q20_dedup.bam',
    output:
        './sMNseq/pilot/samples/{sample}/sorted_{sample}_dedup_unique_read.count'
    log: "./sMNseq/pilot/log/{sample}_count.log"
    resources:
        mem="32G",
        time="04:00:00",
    threads: 3
    shell:
        "samtools view -c {input.bam} > {output} 2>{log}"




####multiBamSummary bins -b **/sorted_*_dedup.bam -o /fast/projects/NB_CircleSeq/work/Micronuclei_sequencing/cov_all_samples.txt --smartLabels -bs 1000000 max/2 --outRawCounts --minMappingQuality 20 --centerReads 