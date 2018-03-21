#import os, sys
import subprocess
configfile: "config.yaml"
(Samples,) = glob_wildcards("../fastq/{sample}_R1_001.fastq.gz")

rule all:
    input:
        #expand("../bam/{sample}_coordsorted_nochr.bam.bai" , sample=Samples),
        #expand("../fastqc/{sample}_R2_001_fastqc.html" , sample=Samples),
        #expand("../fastqc/{sample}_R1_trim_fastqc.html" , sample=Samples),
        expand("../CovMetrics/{sample}_HSmetrics.txt" , sample=Samples),
        #expand("../fastqc/{sample}_R1_trim_fastqc.html", sample=Samples),
        expand("../bam/{sample}_coordsorted_nochr.bam", sample=Samples),

rule fastqc:
    input:
        fq1="../fastq/{sample}_R1_001.fastq.gz",
        fq2="../fastq/{sample}_R2_001.fastq.gz",
    output:
        html1="../fastqc/{sample}_R1_001_fastqc.html",
        html2="../fastqc/{sample}_R2_001_fastqc.html",
        zip1=temp("../fastqc/{sample}_R1_001_fastqc.zip"),
        zip2=temp("../fastqc/{sample}_R2_001_fastqc.zip"),
    threads:config['all']['THREADS']
    params:
        fastqc_dir="../fastqc/"
    conda:
        "envs/fastqc.yaml"
    log:"../logs/fastqc/{sample}.log"
    shell:
        "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"  #Java??

rule SeqPurge_trimming:
    input:
        fq1="../fastq/{sample}_R1_001.fastq.gz",
        fq2="../fastq/{sample}_R2_001.fastq.gz",
    output:
        out1="../trimmed/{sample}_R1_trim.fq.gz",
        out2="../trimmed/{sample}_R2_trim.fq.gz",
    threads:config['all']['THREADS']
    params:
        min_len=config["trim"]["min_len"],
        SeqPurge=config["trim"]["SeqPurge"],
    log:"../logs/trimming/{sample}.log"
    shell:
        "{params.SeqPurge} -min_len 20 -threads {threads} -in1 {input.fq1} -in2 {input.fq2} -out1 {output.out1} -out2 {output.out2}"

rule fastqc_trimmed:
    input:
        fq1="../trimmed/{sample}_R1_trim.fq.gz",
        fq2="../trimmed/{sample}_R2_trim.fq.gz",
    output:
        html1="../fastqc/{sample}_R1_trim_fastqc.html",
        html2="../fastqc/{sample}_R2_trim_fastqc.html",
        zip1=temp("../fastqc/{sample}_R1_trim_fastqc.zip"),
        zip2=temp("../fastqc/{sample}_R2_trim_fastqc.zip"),
    threads:config['all']['THREADS']
    conda:
        "envs/fastqc.yaml"
    params:
        fastqc_dir="../fastqc/"
    log:"../logs/fastqc/{sample}.trim.log"
    shell:
        "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"
        #"fastqc -j java  --threads {threads} {input.fq1} {input.fq2}"  #Java??

def getRGinfo(wildcards):
    fq="../trimmed/"+wildcards.sample+"_R1_trim.fq.gz"
    flowcell=subprocess.check_output("zcat "+fq+" | head -n 1  | awk 'BEGIN {FS=\":\"} {print $3}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()

    lane=subprocess.check_output("zcat "+fq+" | head -n 1 | awk 'BEGIN { FS=\":\"} {print $4}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
    RGinfo="\'@RG\tID:"+lane+"\tSM:"+name+"\tPL:ILLUMINA\tLB:SeqCap\tPU:"+flowcell+"\'"

    return(RGinfo)

rule bwa_mem:
    input:
        fqtrim1="../trimmed/{sample}_R1_trim.fq.gz",
        fqtrim2="../trimmed/{sample}_R2_trim.fq.gz",
    output:
         temp("../bam/{sample}_aligned_reads.bam"),
    params:
        ref= config['all']['REF_CHR'],
        RGinfo=getRGinfo
    threads: config['all']['THREADS']
    conda:
        "envs/bwa-mem.yaml"
    log: "../logs/bwa/{sample}.log"
    shell:
        "bwa mem -M -t {threads} -R {params.RGinfo} {params.ref} {input.fqtrim1} {input.fqtrim2} 2> {log}| samtools view -b -@ {threads}  -> {output} 2>> {log}"

        #"bwa aln -n {params.n} -t {threads} -q {params.q} {params.ref} {input} > {output.sai} 2> {log};"
        #"bwa samse -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -f {output.samse}"
        #" {paramAbra=$JobName"_abra-realigned.bam"s.ref} {output.sai} {input} 2>> {log}"

rule Sambamba_sort:
    input:
        "../bam/{sample}_aligned_reads.bam",
    output:
        sortbam=temp("../bam/{sample}_sorted_reads.bam"),
        sortbai=temp("../bam/{sample}_sorted_reads.bam.bai"),
    params:
        tmpdir=temp("../tmp"),
    conda:
        "envs/sambamba.yaml"
    threads: config['all']['THREADS'],
    log:"../logs/sambamba/{sample}_sort.log",
    shell:
       "sambamba sort -o {output.sortbam} --tmpdir {params.tmpdir} -t {threads} {input} 2> {log} &&"
	   "sambamba index -t {threads} {output.sortbam} {output.sortbai} 2>> {log}"

# rule Samtools_sort:
#     input:
#         "../bam/{sample}_aligned_reads.bam",
#     output:
#         sortbam="../bam/{sample}_sorted_reads.bam",
#         sortbai="../bam/{sample}_sorted_reads.bam.bai",
#     params:
#         tmpdir="../abra/{sample}/tmp",
#     threads: config['all']['THREADS'],
#     log:"../logs/samtools/{sample}_sort.log",
#     shell:
#         "samtools sort -o {output.sortbam} -@ {threads} {input} 2>{log} &&"
# 	    "samtools index {output.sortbam} 2>>{log}"

rule ABRA:
    input:
        sortbam="../bam/{sample}_sorted_reads.bam",
        sortbai="../bam/{sample}_sorted_reads.bam.bai",
    output:
        sv="../abra/{sample}_Abra.sv.txt",
        abra=temp("../bam/{sample}_abra-realigned.bam"),
    params:
        ref= config['all']['REF_CHR'],
        targets=config['all']['targets'],
        tmpdir="../abra/{sample}/",
        abra_prog=config["abra"]["abra_prog"]
    threads: config['all']['THREADS'],
    log:"../logs/abra/{sample}.log",
    shell:
        "java -Xmx8G -jar {params.abra_prog} --in {input.sortbam} --out {output.abra} "
        "--ref {params.ref} --targets {params.targets} --threads {threads} --working {params.tmpdir} --sv {output.sv} &> {log}"

rule mark_duplicates:
    input:
        abra="../bam/{sample}_abra-realigned.bam",
    output:
        qsorted=temp("../bam/{sample}_querysorted.bam"),
        dedup=temp("../bam/{sample}_dedupped.bam"),
        metrics_file="../CovMetrics/{sample}.metrics.txt",
    log: "../logs/mark_duplicates/{sample}.log",
    params:
        tmpdir=temp("../tmp")
    threads: config['all']['THREADS'],
    conda:
        "envs/picard.yaml"
    shell:
        "sambamba sort -o {output.qsorted} --sort-by-name --tmpdir={params.tmpdir} -t {threads} {input.abra} &&"
        "picard MarkDuplicates I={output.qsorted} O={output.dedup} M={output.metrics_file} "
        "ASSUME_SORT_ORDER=queryname &> {log} ;"

rule Sambamba_sort2:
    input:
        dedup="../bam/{sample}_dedupped.bam",
    output:
        coordsorted="../bam/{sample}_coordsorted.bam",
        sortbai="../bam/{sample}_coordsorted.bam.bai",
    params:
        tmpdir=temp("../tmp"),
    threads: config['all']['THREADS'],
    conda:
        "envs/picard.yaml"
    log:"../logs/sambamba/{sample}_sort2.log",
    shell:
        "sambamba sort -o {output.coordsorted} --tmpdir {params.tmpdir} -t {threads} {input.dedup} 2> {log} &&"  #change? samtools
	    "picard BuildBamIndex I={output.coordsorted} 2>> {log}"                      #change? samtools

rule CollectHsMetrics:
    input:
        coordsorted="../bam/{sample}_coordsorted.bam",
    output:
        HSmetrics="../CovMetrics/{sample}_HSmetrics.txt",
        PerTargetCov="../CovMetrics/{sample}_PerTargetCov.txt",
        PerBaseCov="../CovMetrics/{sample}_PerBaseCov.txt",
    params:
        BAIT_INTERVALS=config["all"]["bait_intervals"],
        TARGET_INTERVALS=config["all"]["target_intervals"],
        ref=config["all"]["REF_CHR"],
    conda:
        "envs/picard.yaml"
    threads: config['all']['THREADS'],
    log: "../logs/picard/{sample}.metrics.log",
    shell:
        "picard CollectHsMetrics BAIT_INTERVALS={params.BAIT_INTERVALS} TARGET_INTERVALS={params.TARGET_INTERVALS} INPUT={input.coordsorted} "
        "OUTPUT={output.HSmetrics} METRIC_ACCUMULATION_LEVEL=ALL_READS PER_TARGET_COVERAGE={output.PerTargetCov} "
        "PER_BASE_COVERAGE={output.PerBaseCov} REFERENCE_SEQUENCE={params.ref} 2> {log}"

rule remove_chr_prefix :
    input:
        coordsorted="../bam/{sample}_coordsorted.bam",
    output:
        nochrsam=temp("../bam/{sample}_coordsorted_nochr.sam"),
        nochrbam="../bam/{sample}_coordsorted_nochr.bam",
        nochrbai="../bam/{sample}_coordsorted_nochr.bam.bai",
    params:
    conda:
        "envs/samtools.yaml"
    threads: config['all']['THREADS'],
    log: "../logs/samtools/{sample}-remove_prefix.txt"
    shell:
        "samtools view -@ {threads} -h {input.coordsorted} | sed -e 's/chr//g' > {output.nochrsam} && "
        "samtools view -@ {threads} -Sb  {output.nochrsam} > {output.nochrbam} && "
        "samtools index {output.nochrbam} "

rule BaseRecalibration:
        input:
                coordsorted="../bam/{sample}_coordsorted.bam",
        output:
                recal_table="../baseRecal/{sample}_recal_data.table",
                bam="../bam/{sample}_recal.bam",
        params:
                ref=config["all"]["REF_CHR"],
                targets=config['all']['targets'],
                known_indels=config['baseRecalibration']['known_indels'],
                known_snps=config['baseRecalibration']['known_snps'],
        threads: config['all']['THREADS']
        conda:
                "envs/gatk4.yaml"
        log: "../logs/baseRecalibration/{sample}_baseRecalibration.txt"
        shell:
                "gatk-launch BaseRecalibrator -R {params.ref} -I {input.coordsorted} -L {params.targets} "
                "-known-sites {params.known_indels} -known-sites {params.known_snps} --output {output.recal_table} &> {log} &&"
                "gatk-launch ApplyBQSR -R {params.ref} -I {input.coordsorted} --bqsr-recal-file {output.recal_table} -O {output.bam} &>> {log}"

# rule :
#     input:
#
#     output:
#
#     params:
#
#     threads: config['all']['THREADS']
#     log:
#     shell:
