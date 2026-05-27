#Snakemake workflow for trimming and mapping ChIPseq data. 

#Current genome: hg38 (canonical_hg38)

##snakemake --snakefile run_bwa_mapping_hg38.smk --profile=cubi-v1 -j 18 -k -p --latency-wait 120 --rerun-incomplete -n


SAMPLES = ["A4714_1_DM_HU_Ac_S4_L001","A4714_2_DM_UT_Ac_S5_L001","A4714_10_DM_UT_INPUT_S1_L001"]


print(SAMPLES)


#define controls
sample_inputcontrols = {
  "A4714_1_DM_HU_Ac_S4_L001": ["A4714_10_DM_UT_INPUT_S1_L001"],
  "A4714_2_DM_UT_Ac_S5_L001": ["A4714_10_DM_UT_INPUT_S1_L001"],
  "A4714_10_DM_UT_INPUT_S1_L001": ["A4714_10_DM_UT_INPUT_S1_L001"]
}

INPUT_DIR =""
OUT = ""
LOG = ""

rule all:
  input:
    expand(OUT + '/{sample}/{sample}_R1_001_trimmed.fq.gz', sample=SAMPLES),
    expand(LOG + '/{sample}_fastqc.log', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted.bam', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}.flagstat', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered.bam', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup.bam', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup.stats', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam.bai', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.flagstat', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final_10bp_10Mb.bw', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.Input_10bp_10Mb_subtracted.bw', sample=SAMPLES),
    expand(OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.phantompeakqualtools.txt', sample=SAMPLES)



rule trim:
  threads: 8
  resources:
    mem='64G',
    time='32:00:00',
  input:
    read1=INPUT_DIR + '/{sample}_R1_001.fastq.gz'
  params:
    out_dir=OUT + '/{sample}/'
  output:
    OUT + '/{sample}/{sample}_R1_001_trimmed.fq.gz'
  log: LOG + '/{sample}_trim.log'
  shell:
    'trim_galore --nextseq 20 --phred33 --illumina -j 8 -o {params.out_dir} {input.read1} 2>{log}'

rule fastqc:
  threads: 8
  resources:
    mem='64G',
    time='32:00:00'
  input:
    read1=OUT + '/{sample}/{sample}_R1_001_trimmed.fq.gz'
  params:
    out_dir=OUT + '/{sample}/',
    fastqs=OUT + '/{sample}/{sample}*fq'
  log:
    LOG + '/{sample}_fastqc.log'
  shell:
    'fastqc -t 6 -o {params.out_dir} {params.fastqs} 2>{log}> '



rule bwa:
  input:
    fastq=OUT + '/{sample}/{sample}_R1_001_trimmed.fq.gz'
  params:
    group='"@RG\\tID:Histone_ChIP_ecDNA_MN\\tPL:ILLUMINA\\tLB:{sample}\\tSM:{sample}"',
    genome="./ref/GRCh38.primary_assembly.canonical.EBV.fa"
  output:
    sam=temp(OUT + '/{sample}/canonical_hg38/{sample}.sam')
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log:
    LOG + '/{sample}_bwa.log'
  shell:
    'bwa mem -q -t {threads} -R {params.group} {params.genome} {input.fastq} > {output.sam} 2>{log}'


rule coordinate_sort1:
  input:
    sam=OUT + '/{sample}/canonical_hg38/{sample}.sam'
  output:
    sorted_bam=OUT + '/{sample}/canonical_hg38/{sample}_sorted.bam'
  resources:
    mem="32G",
    time="12:00:00",
  threads: 8
  log: LOG + '/{sample}_coord_sort1.log'
  shell:
    'samtools sort -@ {threads} -o {output.sorted_bam} {input.sam} 2>{log}'



rule flagstat1:
  input:
    bam_sorted=OUT + '/{sample}/canonical_hg38/{sample}_sorted.bam'
  output:
    flagstat=OUT + '/{sample}/canonical_hg38/{sample}.flagstat'
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_flagstat.log'
  shell:
    'samtools flagstat {input.bam_sorted} > {output.flagstat} 2>{log}'
 

rule read_filter:
  input:
    bam_sorted=OUT + '/{sample}/canonical_hg38/{sample}_sorted.bam'
  output:
    bam_sorted_filtered=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered.bam'
  params:
    mapq = 20
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_read_filter.log'
  shell:
    'samtools view -F 1804 -q {params.mapq} -b {input.bam_sorted} -o {output.bam_sorted_filtered} 2>{log}'



rule mark_duplicates:
  input:
    bam_sorted_filtered=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered.bam'
  output:
    bam_sorted_filtered_mkdup=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup.bam',
    qc=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup.stats'
  params: 
    mapq = 20
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_markdup.log'
  shell:
    'picard MarkDuplicates INPUT={input.bam_sorted_filtered} OUTPUT={output.bam_sorted_filtered_mkdup} METRICS_FILE={output.qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false 2>{log}'



rule rm_dups:
  input:
    bam_sorted_filtered_mkdup=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup.bam'
  output:
    bam_sorted_filtered_mkdup_rm=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam'
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_read_filter2final.log'
  shell:
    'samtools view -F 1796 -b {input.bam_sorted_filtered_mkdup} -o {output.bam_sorted_filtered_mkdup_rm} 2>{log}'

rule index:
  input:
    bam=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam'
  output:
    bai=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam.bai'
  resources:
    mem="61G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_index.log'
  shell:
    'samtools index {input.bam} 2>{log}'


rule flagstat2:
  input:
    bam_sorted=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam'
  output:
    flagstat=OUT + '/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.flagstat'
  resources:
    mem="64G",
    time="12:00:00"
  threads: 10
  log: LOG + '/{sample}_flagstat2.log'
  shell:
    'samtools flagstat {input.bam_sorted} > {output.flagstat} 2>{log}'



rule bigWig:
  input:
    OUT + "/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam"
  output:
    OUT +"/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final_10bp_10Mb.bw"
  params:
    blacklist = "./reference/blacklisted_regions_hg38.bed"
  threads: 8
  resources:
    mem='32G',
    time='12:00:00',
    partition='medium'
  shell:
    "bamCoverage --bam {input} -o {output} --binSize 10 --ignoreForNormalization chrX chrM NC_007605.1 --scaleFactor 10 --effectiveGenomeSize 2805636231 --exactScaling --extendReads 200 --blackListFileName {params.blacklist} --normalizeUsing CPM"

rule bigWig_FC:
  input:
    treatment = OUT + "/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final_10bp_10Mb.bw",
    control = lambda wildcards: OUT + "/" + sample_inputcontrols[wildcards.sample][0] + "/canonical_hg38/" + sample_inputcontrols[wildcards.sample][0] + "_sorted_filtered_mkdup_final_10bp_10Mb.bw"
  output:
    OUT + "/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.Input_10bp_10Mb_subtracted.bw"
  threads: 8
  resources:
    mem='96G',
    time='12:00:00',
    partition='medium'
  shell:
    "bigwigCompare --bigwig1 {input.treatment} --bigwig2 {input.control} --operation subtract --skipZeroOverZero --skipNAs -p 8 --binSize 10 --pseudocount 0.01 --outFileName {output}"

rule phantompeakqualtools:
  input:
    OUT + "/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.bam"
  output:
    OUT + "/{sample}/canonical_hg38/{sample}_sorted_filtered_mkdup_final.phantompeakqualtools.txt"
  threads: 8
  resources:
    mem='32G',
    time='96:00:00',
    partition='medium'
  shell:
    "run_spp.R -rf -c={input} -savp -out={output}"


#bigwigCompare --bigwig1 {input.treatment} --bigwig2 {input.control} --operation subtract --skipZeroOverZero --skipNAs -p 8 --binSize 10 --pseudocount 0.01 --outFileName {output}
