#import os, sys
import subprocess
import pandas as pd
configfile: "config.yaml"
(Samples,) = glob_wildcards("../fastq/{sample}_R1_001.fastq.gz")


all_Samples=pd.read_csv('matching_samples.tsv', sep='\t')

Tumor=list(all_Samples['Tumor'])
Normal=list(all_Samples['Normal'])

def getnames(samplelist):
    SAMPLES=dict()
    for sample in samplelist:
        matching = [name for name in Samples if sample+"_" in name]  # get the name including lanenumber
        SAMPLES[sample]=matching[0]
		#sample = re.match('[a-zA-Z0-9\-]*', name).group(0)
        #fastqfile=SRdir+"fastq/"+name+".fastq.gz"
        #SAMPLES[sample]=fastqfile
    return(SAMPLES)

Tumor_samples=getnames(Tumor)  # dictionary containing samplename as key and samplename with lanenumber as value
Normal_samples=getnames(Normal) # dictionary containing samplename as key and samplename with lanenumber as value
pairs=dict(zip(Tumor, Normal)) # dictionary containing Tumorname as key and normalname as value

rule all:
    input:
        expand("../Mutect2/{tumor}.somatic.vcf.gz", tumor=Tumor_samples.keys()),
        expand("../Mutect2/{tumor}.somatic_filtered.vcf.gz", tumor=Tumor_samples.keys()),
        #expand("../bam/{sample}_coordsorted_nochr.bam.bai" , sample=Samples),
        expand("../fastqc/{sample}_R1_001_fastqc.html" , sample=Samples),
        expand("../fastqc/{sample}_R1_trim_fastqc.html" , sample=Samples),
        expand("../CovMetrics/{sample}_HSmetrics.txt" , sample=Samples),
        expand("../vcf/annotated/{tumor}.annotated_effect.vcf", tumor=Tumor_samples.keys()),
        expand("../vcf/filtered/{tumor}.LOW.csv", tumor=Tumor_samples.keys())
        #expand("../fastqc/{sample}_R1_trim_fastqc.html", sample=Samples),
        #expand("../bam/{sample}_coordsorted_nochr.bam", sample=Samples),

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

rule Cutadapt_trimming:
    input:
        fq1="../fastq/{sample}_R1_001.fastq.gz",
        fq2="../fastq/{sample}_R2_001.fastq.gz",
    output:
        out1="../trimmed/{sample}_R1_trim.fq.gz",
        out2="../trimmed/{sample}_R2_trim.fq.gz",
    threads: config['all']['THREADS']
    params:
        adapter3="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",   #standard ILLUMINA but also IDT
        adapter5="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT",  #standard ILLUMINA but also IDT
        min_len=config["trim"]["min_len"],
    log: "../logs/trimming/{sample}.log",
    conda:
        "envs/cutadapt.yaml"
    shell:
        "cutadapt -a {params.adapter3} -A {params.adapter5} -j {threads} --minimum-length {params.min_len} "
        " -o {output.out1} -p {output.out2} {input.fq1} {input.fq2}  > {log}"

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
        ref= config['all']['REF'],
        RGinfo=getRGinfo
    threads: config['all']['THREADS']
    conda:
        "envs/bwa-mem.yaml"
    log: "../logs/bwa/{sample}.log"
    shell:
        "bwa mem -M -t {threads} -R {params.RGinfo} {params.ref} {input.fqtrim1} {input.fqtrim2} 2> {log} | samtools view -b -@ {threads}  -> {output} 2>> {log}"

        #"bwa aln -n {params.n} -t {threads} -q {params.q} {params.ref} {input} > {output.sai} 2> {log};"
        #"bwa samse -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -f {output.samse}"
        #" {paramAbra=$JobName"_abra-realigned.bam"s.ref} {output.sai} {input} 2>> {log}"

rule Sambamba_sort:
    input:
        "../bam/{sample}_aligned_reads.bam",
    output:
        qsorted=temp("../bam/{sample}_querysorted.bam"),
        #sortbai=temp("../bam/{sample}_querysorted.bam.bai"),
    params:
        tmpdir=temp("../tmp"),
    conda:
        "envs/sambamba.yaml"
    threads: config['all']['THREADS'],
    log: "../logs/sambamba/{sample}_sort.log",
    shell:
       "sambamba sort -o {output.qsorted} --sort-by-name --tmpdir={params.tmpdir} -t {threads} {input} 2> {log} "


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

# rule ABRA:
#     input:
#         sortbam="../bam/{sample}_sorted_reads.bam",
#         sortbai="../bam/{sample}_sorted_reads.bam.bai",
#     output:
#         sv="../abra/{sample}_Abra.sv.txt",
#         abra=temp("../bam/{sample}_abra-realigned.bam"),
#     params:
#         ref= config['all']['REF'],
#         targets=config['all']['targets'],
#         tmpdir="../abra/{sample}/",
#         abra_prog=config["abra"]["abra_prog"]
#     threads: config['all']['THREADS'],
#     log:"../logs/abra/{sample}.log",
#     shell:
#         "java -Xmx8G -jar {params.abra_prog} --in {input.sortbam} --out {output.abra} "
#         "--ref {params.ref} --targets {params.targets} --threads {threads} --working {params.tmpdir} --sv {output.sv} &> {log}"

rule mark_duplicates:
    input:
        qsorted="../bam/{sample}_querysorted.bam",
    output:
        dedup=temp("../bam/{sample}_dedupped.bam"),
        metrics_file="../CovMetrics/{sample}.metrics.txt",
    log: "../logs/mark_duplicates/{sample}.log",
    params:
        tmpdir=temp("../tmp")
    threads: config['all']['THREADS'],
    conda:
        "envs/picard.yaml"
    shell:
        "picard -Xms1g -Xmx8g MarkDuplicates I={input.qsorted} O={output.dedup} M={output.metrics_file} "
        "ASSUME_SORT_ORDER=queryname TMP_DIR={params.tmpdir} &> {log} ;"

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
        "sambamba sort -o {output.coordsorted} --tmpdir {params.tmpdir} -t {threads} {input.dedup} 2> {log} && "
        "sambamba index -t {threads} {output.coordsorted} {output.sortbai} 2>> {log}"
        #"picard BuildBamIndex I={output.coordsorted} 2>> {log}"                      #change? samtools


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
        ref=config["all"]["REF"],
    conda:
        "envs/picard.yaml"
    threads: config['all']['THREADS'],
    log: "../logs/picard/{sample}.metrics.log",
    shell:
        "picard -Xmx4g CollectHsMetrics "
        "BAIT_INTERVALS={params.BAIT_INTERVALS} "
        "TARGET_INTERVALS={params.TARGET_INTERVALS} "
        "INPUT={input.coordsorted} "
        "OUTPUT={output.HSmetrics} "
        "METRIC_ACCUMULATION_LEVEL=ALL_READS "
        "PER_TARGET_COVERAGE={output.PerTargetCov} "
        "PER_BASE_COVERAGE={output.PerBaseCov} "
        "REFERENCE_SEQUENCE={params.ref} 2> {log}"

rule Mutect2:
	input:
		tumor=lambda wildcards:"../bam/"+Tumor_samples[wildcards.tumor]+"_coordsorted.bam",
		normal=lambda wildcards:"../bam/"+Normal_samples[pairs[wildcards.tumor]]+"_coordsorted.bam"
	output:
		vcf="../Mutect2/{tumor}.somatic.vcf.gz",
		bamout="../Mutect2/{tumor}-N_m2.bam"
	params:
		ref=config["all"]["REF"],
		snps=config["all"]["snps"],
		normalname=lambda wildcards:pairs[wildcards.tumor]
	threads: config['all']['THREADS']
	conda:
		"envs/gatk4.yaml"
	log: "../logs/Mutect2/{tumor}_mutect2.txt"
	shell:
		"gatk Mutect2 -R {params.ref} -I {input.tumor} -tumor {wildcards.tumor} -I {input.normal} -normal"
		" {params.normalname} --germline-resource {params.snps} -O {output.vcf} -bamout {output.bamout} "
		"--native-pair-hmm-threads {threads} &> {log}"

#Specify the case sample for somatic calling with two parameters. Provide the BAM with -I and the sample's read group sample name (the SM field value) with -tumor. To look up the read group SM field use GetSampleName. Alternatively, use samtools view -H tumor.bam | grep '@RG'.

rule filtermutect2:
	input:
		vcf="../Mutect2/{tumor}.somatic.vcf.gz",
	output:
		filt="../Mutect2/{tumor}.somatic_filtered.vcf.gz"
	conda:
		"envs/gatk4.yaml"
	log: "../logs/Mutect2/{tumor}_filtermutect2.txt"
	shell:
		"gatk FilterMutectCalls -V {input.vcf} -O {output.filt} &> {log} "


rule SnpSift:
	input:
		filt="../Mutect2/{tumor}.somatic_filtered.vcf.gz"
	output:
		annotated="../vcf/annotated/{tumor}.annotated.vcf"
	conda:
		"envs/SnpSift.yaml"
	params:
		Java_mem=config["all"]["Java_mem"],
		dbsnp=config["all"]["dbsnp"],
		clinvar=config["all"]["clinvar"],
		Cosmic=config["all"]["Cosmic"]
	log:
		"../logs/snpEff/{tumor}_annotation.txt"
	shell:
		"SnpSift annotate {params.Java_mem} {params.dbsnp} {input.filt} 2> {log} | "
		"SnpSift annotate {params.Java_mem} {params.clinvar} -  2>> {log} | "
		"SnpSift annotate {params.Java_mem} {params.Cosmic} - > {output.annotated} 2>> {log} "

rule SnpEff:
	input:
		annotated="../vcf/annotated/{tumor}.annotated.vcf"
	output:
		effect="../vcf/annotated/{tumor}.annotated_effect.vcf",
		snpEff_stats="../vcf/annotated/stats/{tumor}_stats.html"
	conda:
		"envs/snpeff.yaml"
	params:
		Java_mem=config["all"]["Java_mem"],
		targets=config["all"]["targets"]
	log:
		"../logs/snpEff/{tumor}_snpEff.txt"
	shell:
		"snpEff {params.Java_mem} eff -filterInterval {params.targets} -v -canon -strict "
		"-stats {output.snpEff_stats} hg19 {input.annotated} > {output.effect} 2> {log} "


rule SnpSift_filter:
	input:
		effect="../vcf/annotated/{tumor}.annotated_effect.vcf",
	output:
		low="../vcf/filtered/{tumor}.LOW.vcf",
		high="../vcf/filtered/{tumor}.HIGH.vcf"
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift filter -f {input.effect} "(FILTER='PASS')&((ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE'))"> {output.high} && \
		SnpSift filter -f {input.effect} "(FILTER='PASS')" | SnpSift filter -n "(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > {output.low}
		"""

rule SnpSift_csv:
	input:
		high="../vcf/filtered/{tumor}.HIGH.vcf",
		low="../vcf/filtered/{tumor}.LOW.vcf"
	output:
		high="../vcf/filtered/{tumor}.HIGH.csv",
		low="../vcf/filtered/{tumor}.LOW.csv"
	params:
		fields='CHROM POS "ANN[0].GENE" REF ALT DP "GEN[0].AD" "GEN[1].AD" "GEN[0].AF" "GEN[1].AF" "ANN[0].IMPACT" \
		"ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" ID "COMMON" "RS" POP_AF "CAF"  "LOF" "NMD" "MUT" \
		"CLNSIG" "ORIGIN" "SNP" "COSM.ID" "FATHMM" "MUT.ST"'
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift extractFields -e "." {input.high} {params.fields} > {output.high} && \
		SnpSift extractFields -e "." {input.low} {params.fields} > {output.low}
		"""
#GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB
### Create csv files - extract fields of interest
#java -jar $snpEff/SnpSift.jar extractFields -e "." $vcf_out \
#CHROM POS "ANN[0].GENE" REF ALT DP AF \
#"ANN[0].IMPACT" "ANN[0].EFFECT" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
#ID "COMMON" "RS" "CAF" "LOF" "NMD" "MUT" \
#"CLNSIG" "ORIGIN" "SNP" "AF_EXAC" "AF_TGP" "gnomAD_AF" \
#"COSM.ID" "FATHMM" "MUT.ST" "PON_BLACKLIST"> $vcf_csv

#	java -jar $SnpSift filter -f $somatic_vcf \
#	"(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > $HI_somatic_vcf

#	java -jar $SnpSift filter -f $somatic_vcf -n \
#	"(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > $LOW_somatic_vcf



# rule remove_chr_prefix :
#     input:
#         coordsorted="../bam/{sample}_coordsorted.bam",
#     output:
#         nochrsam=temp("../bam/{sample}_coordsorted_nochr.sam"),
#         nochrbam="../bam/{sample}_coordsorted_nochr.bam",
#         nochrbai="../bam/{sample}_coordsorted_nochr.bam.bai",
#     params:
#     conda:
#         "envs/samtools.yaml"
#     threads: config['all']['THREADS'],
#     log: "../logs/samtools/{sample}-remove_prefix.txt"
#     shell:
#         "samtools view -@ {threads} -h {input.coordsorted} | sed -e 's/chr//g' > {output.nochrsam} && "
#         "samtools view -@ {threads} -Sb  {output.nochrsam} > {output.nochrbam} && "
#         "samtools index {output.nochrbam} "

# rule BaseRecalibration:
#         input:
#                 coordsorted="../bam/{sample}_coordsorted.bam",
#         output:
#                 recal_table="../baseRecal/{sample}_recal_data.table",
#                 bam="../bam/{sample}_recal.bam",
#         params:
#                 ref=config["all"]["REF"],
#                 targets=config['all']['targets'],
#                 known_indels=config['baseRecalibration']['known_indels'],
#                 known_snps=config['baseRecalibration']['known_snps'],
#         threads: config['all']['THREADS']
#         conda:
#                 "envs/gatk4.yaml"
#         log: "../logs/baseRecalibration/{sample}_baseRecalibration.txt"
#         shell:
#                 "gatk-launch BaseRecalibrator -R {params.ref} -I {input.coordsorted} -L {params.targets} "
#                 "-known-sites {params.known_indels} -known-sites {params.known_snps} --output {output.recal_table} &> {log} &&"
#                 "gatk-launch ApplyBQSR -R {params.ref} -I {input.coordsorted} --bqsr-recal-file {output.recal_table} -O {output.bam} &>> {log}"


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
