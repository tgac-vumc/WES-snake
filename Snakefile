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

#Tumor_samples=getnames(Tumor)  # dictionary containing samplename as key and samplename with lanenumber as value
#Normal_samples=getnames(Normal) # dictionary containing samplename as key and samplename with lanenumber as value

# I merged samples sequenced first time and second time after deduplication by hand (merge_bam.sh) for the merged samples I don't use the lanenumber, so I don't need getnames, because I don't want to change the whole script I keep using a dict format.
Tumor_samples=dict(zip(Tumor, Tumor))
Normal_samples=dict(zip(Normal, Normal))
pairs=dict(zip(Tumor, Normal)) # dictionary containing Tumorname as key and normalname as value

rule all:
    input:
        expand("../Mutect2/{tumor}.somatic.vcf.gz", tumor=Tumor_samples.keys()),
        expand("../Mutect2/{tumor}.somatic_filtered.vcf.gz", tumor=Tumor_samples.keys()),
        #expand("../bam/{sample}_coordsorted_nochr.bam.bai" , sample=Samples),
        #expand("../fastqc/{sample}_R1_001_fastqc.html" , sample=Samples),
        #expand("../fastqc/{sample}_R1_trim_fastqc.html" , sample=Samples),
        expand("../CovMetrics/{sample}_HSmetrics.txt" , sample=Tumor_samples.keys()),
        expand("../CovMetrics/{sample}_HSmetrics.txt" , sample=Normal_samples.keys()),
        expand("../vcf/annotated/{tumor}.annotated_effect.vcf", tumor=Tumor_samples.keys()),
        expand("../vcf/filtered/{tumor}.all_evidenced.csv", tumor=Tumor_samples.keys()),
        expand("../vcf/filtered/{tumor}.all_evidenced.vcf", tumor=Tumor_samples.keys()),
		#"../MutationalPatterns/mutationspectrum_96.png"
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
        TARGET_INTERVALS=config["all"]["bait_intervals"],
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
		normalname=lambda wildcards:pairs[wildcards.tumor],
		tmpdir=temp("../tmp"),
	threads: config['all']['THREADS']
	conda:
		"envs/gatk4.yaml"
	log: "../logs/Mutect2/{tumor}_mutect2.txt"
	shell:
		"gatk Mutect2 -R {params.ref} -I {input.tumor} -tumor {wildcards.tumor} -I {input.normal} -normal"
		" {params.normalname} --germline-resource {params.snps} -O {output.vcf} -bamout {output.bamout} "
		"--native-pair-hmm-threads {threads} --TMP_DIR {params.tmpdir} &> {log}"

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

#TODO make params of SnpSift_filters
rule SnpSift_filter:
	input:
		effect="../vcf/annotated/{tumor}.annotated_effect.vcf",
	output:
		highimpact="../vcf/filtered/{tumor}.high_impact.vcf",
		lowimpact="../vcf/filtered/{tumor}.low_impact.vcf",
		all_evidenced="../vcf/filtered/{tumor}.all_evidenced.vcf"
	params:
		minimal_Alternative_depth=4
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift filter -f {input.effect} "(GEN[1].AD[1] >= 4) & (GEN[0].AD[0] + GEN[0].AD[1] > 10) & (GEN[1].AD[0] + GEN[1].AD[1] > 10) & (FILTER='PASS') &  (GEN[1].F1R2[1] >= 1) & (GEN[1].F2R1[1] >= 1) & GEN[1].AF >= 0.05" > {output.all_evidenced} && \
		SnpSift filter -f {output.all_evidenced} -n "(ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE')" > {output.lowimpact} && \
		SnpSift filter -f {output.all_evidenced} "((ANN[0].IMPACT='HIGH') | (ANN[0].IMPACT='MODERATE'))" > {output.highimpact}
		"""

rule SnpSift_csv:
	input:
		highimpact="../vcf/filtered/{tumor}.high_impact.vcf",
		lowimpact="../vcf/filtered/{tumor}.low_impact.vcf",
		all_evidenced="../vcf/filtered/{tumor}.all_evidenced.vcf"
	output:
		highimpact="../vcf/filtered/{tumor}.high_impact.csv",
		lowimpact="../vcf/filtered/{tumor}.low_impact.csv",
		all_evidenced="../vcf/filtered/{tumor}.all_evidenced.csv"
	params:
		fields='CHROM POS "ANN[0].GENE" REF ALT DP "GEN[0].AD" "GEN[1].AD" "GEN[0].AF" "GEN[1].AF" "ANN[0].IMPACT" \
		"ANN[0].EFFECT" FILTER "GEN[1].F1R2" "GEN[1].F2R1" "GEN[1].MBQ" "GEN[1].MMQ" "ANN[0].HGVS_C" \
		"ANN[0].HGVS_P" ID "COMMON" "RS" POP_AF "CAF"  "LOF" "NMD" "MUT" \
		"CLNSIG" "ORIGIN" "SNP" "COSM.ID" "FATHMM" "MUT.ST"'
	conda:
		"envs/SnpSift.yaml"
	shell:
		"""
		SnpSift extractFields -e "." {input.highimpact} {params.fields} > {output.highimpact} && \
		SnpSift extractFields -e "." {input.lowimpact} {params.fields} > {output.lowimpact} && \
		SnpSift extractFields -e "." {input.all_evidenced} {params.fields} > {output.all_evidenced}
		"""

# rule MutationalPatterns:
# 	input:
# 		all_evidenced="../vcf/filtered/{tumor}.all_evidenced.vcf",
# 		script="scripts/Mutational_patterns.R"
# 	output:
# 		figure="../MutationalPatterns/mutationspectrum_96.png"
# 	conda:
# 		"envs/R.yaml"
# 	script:
# 		"{input.script}"
