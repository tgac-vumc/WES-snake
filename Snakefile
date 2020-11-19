#import os, sys
import subprocess
import pandas as pd
import glob
import os.path as path
configfile: "config.yaml"

## configurations ###########
PATH_FASTQ = config['path']['fastq']
PATH_QC = config['path']['qc']
PATH_BAM = config['path']['bam']
PATH_LOG = config['path']['log']
PATH_TEMP = config['path']['temp']
PATH_PIPELINE = srcdir('.')

PLATFORM = config['platform']['SRorPE'] # SE or PE
PREFIX = config['platform']['prefix'] # prefix before "fastq.gz"
#################################


## obtain sample list ###########
all_Samples=pd.read_csv('test_samples.tsv', sep='\t')

Tumor=list(all_Samples['Tumor']) # these are IDs, and use as the key to fetch files.
Normal=list(all_Samples['Normal'])


def getnames(samplelist, platform, pathdata):
    Files = []
    RNAIDs = []
    SAMPLES = dict()
    for sample in samplelist:
        if platform in ['SR', 'sr']:
            for prefix in PREFIX:
                SAMPLES[sample] = glob.glob(path.join(pathdata, sample+'*'+prefix+'.fastq.gz'))
                SAMPLES[sample].sort()
        if platform in ['PE', 'pe']:
            SAMPLES[sample] = dict()
            SAMPLES[sample]['R1'] = glob.glob(path.join(pathdata, sample+'*'+PREFIX[0]+'.fastq.gz'))
            SAMPLES[sample]['R1'].sort()
            SAMPLES[sample]['R2'] = glob.glob(path.join(pathdata, sample+'*'+PREFIX[1]+'.fastq.gz'))
            SAMPLES[sample]['R2'].sort()
    return(SAMPLES)

Tumor_samples = getnames(Tumor, PLATFORM, PATH_FASTQ)
Normal_samples = getnames(Normal, PLATFORM, PATH_FASTQ)
AllFiles = Tumor_samples
AllFiles.update(Normal_samples)

print(Tumor_samples)
print(Normal_samples)

#Tumor_samples=getnames(Tumor)  # dictionary containing samplename as key and samplename with lanenumber as value
#Normal_samples=getnames(Normal) # dictionary containing samplename as key and samplename with lanenumber as value
# TODO this can be made easier but need alterations in the whole script so this is the quick solution. 
#Tumor_samples=dict(zip(Tumor, Tumor))
#Normal_samples=dict(zip(Normal, Normal)

pairs=dict(zip(Tumor, Normal)) # dictionary containing Tumorname as key and normalname as value


#####################################


rule all:
    input:
#        expand(path.join(PATH_COV, "{sample}_HSmetrics.txt"), sample=Tumor_samples.keys()),
#        expand(path.join(PATH_COV, "{sample}_HSmetrics.txt"), sample=Normal_samples.keys()),
#        path.join(PATH_COV, 'Metrics_allsample.txt')
#       expand(path.join(PATH_QC, 'Combined/{sample}_Combined_metrics.txt'), sample=AllFiles.keys()),
        path.join(PATH_QC, 'Combined_Metrics.txt'),
        expand(path.join(PATH_QC, 'multiqc{trim}.html'), trim = ['', '_trim']) 
        #expand("../vcf/filtered/{tumor}.all_evidenced.vcf", tumor=Tumor_samples.keys())
 

def getRGinfo(wildcards):
    fq=path.join(PATH_FASTQ, "trimmed", wildcards.sample+"_R1_trim.fq.gz")
    flowcell=subprocess.check_output("zcat "+fq+" | head -n 1  | awk 'BEGIN {FS=\":\"} {print $3}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()

    lane=subprocess.check_output("zcat "+fq+" | head -n 1 | awk 'BEGIN { FS=\":\"} {print $4}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
    RGinfo="\'@RG\tID:"+lane+"\tSM:"+name+"\tPL:ILLUMINA\tLB:SeqCap\tPU:"+flowcell+"\'"

    return(RGinfo)


if PLATFORM in ['PE', 'pe']: # I only implement PE for now

    rule merge_fastq:
        input:
            lambda wildcards: AllFiles[wildcards.sample][wildcards.R1R2]
        output:
            path.join(PATH_FASTQ, '{sample}.fastq.merged{R1R2}.gz')
        shell:
            "cat {input} > {output}"


    rule fastqc:
        input:
            fq1=path.join(PATH_FASTQ, "{sample}.fastq.mergedR1.gz"),
            fq2=path.join(PATH_FASTQ, "{sample}.fastq.mergedR2.gz")
        output:
            html1=path.join(PATH_QC, "fastqc", "{sample}.fastq.mergedR1_fastqc.html"),
            html2=path.join(PATH_QC, "fastqc", "{sample}.fastq.mergedR2_fastqc.html"),
            zip1=temp(path.join(PATH_QC, "fastqc", "{sample}.fastq.mergedR1_fastqc.zip")),
            zip2=temp(path.join(PATH_QC, "fastqc", "{sample}.fastq.mergedR2_fastqc.zip"))
        threads:config['all']['THREADS']
        params:
            fastqc_dir=path.join(PATH_QC, "fastqc")
        conda:
            "envs/fastqc.yaml"
        log: path.join(PATH_LOG, "fastqc","{sample}.log")
        shell:
            "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"  #Java??

    rule Cutadapt_trimming:
        input:
            fq1=path.join(PATH_FASTQ, "{sample}.fastq.mergedR1.gz"),
            fq2=path.join(PATH_FASTQ, "{sample}.fastq.mergedR2.gz")
        output:
            out1=path.join(PATH_FASTQ, "trimmed","{sample}_R1_trim.fq.gz"),
            out2=path.join(PATH_FASTQ, "trimmed","{sample}_R2_trim.fq.gz")
        threads: config['all']['THREADS']
        params:
            adapter3=config['cutadapt']['adapter'][0],   #standard ILLUMINA but also IDT
            adapter5=config['cutadapt']['adapter'][1],  #standard ILLUMINA but also IDT
            min_len=config["trim"]["min_len"],
        log: path.join(PATH_LOG, "trimming", "{sample}.log"),
        conda:
            "envs/cutadapt.yaml"
        shell:
            "cutadapt -a {params.adapter3} -A {params.adapter5} -j {threads} --minimum-length {params.min_len} "
            " -o {output.out1} -p {output.out2} {input.fq1} {input.fq2} > {log}"

    rule fastqc_trimmed:
        input:
            fq1=path.join(PATH_FASTQ, "trimmed", "{sample}_R1_trim.fq.gz"),
            fq2=path.join(PATH_FASTQ, "trimmed", "{sample}_R2_trim.fq.gz"),
        output:
            html1=path.join(PATH_QC,"fastqc","{sample}_R1_trim_fastqc.html"),
            html2=path.join(PATH_QC,"fastqc","{sample}_R2_trim_fastqc.html"),
            zip1=temp(path.join(PATH_QC, "fastqc", "{sample}_R1_trim_fastqc.zip")),
            zip2=temp(path.join(PATH_QC, "fastqc", "{sample}_R2_trim_fastqc.zip")),
        threads:config['all']['THREADS']
        conda:
            "envs/fastqc.yaml"
        params:
            fastqc_dir=path.join(PATH_QC, "fastqc")
        log: path.join(PATH_LOG, "fastqc","{sample}.trim.log")
        shell:
            "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"
        #"fastqc -j java  --threads {threads} {input.fq1} {input.fq2}"  #Java??

    rule multiqc_raw:
        input:
            expand(path.join(PATH_QC, "fastqc", "{sample}{pair}_fastqc.zip"),
                            pair=['.fastq.mergedR1','.fastq.mergedR2'], sample=AllFiles.keys())
        output:
            path.join(PATH_QC,'multiqc.html') # empty wildcard allowed
        log:
            path.join(PATH_LOG,'multiqc.log')
        wrapper:
            '0.35.0/bio/multiqc'


    rule multiqc_trim:
        input:
            expand(path.join(PATH_QC, "fastqc", "{sample}{pair}_trim_fastqc.zip"),
                            pair=['_R1','_R2'], sample=AllFiles.keys())
        output:
            path.join(PATH_QC,'multiqc_trim.html') # empty wildcard allowed
        log:
            path.join(PATH_LOG,'multiqc_trim.log')
        wrapper:
            '0.35.0/bio/multiqc'

    rule bwa_mem:
        input:
            fqtrim1=path.join(PATH_FASTQ, "trimmed", "{sample}_R1_trim.fq.gz"),
            fqtrim2=path.join(PATH_FASTQ, "trimmed", "{sample}_R2_trim.fq.gz"),
        output:
             temp(path.join(PATH_BAM, "{sample}_aligned_reads.bam")),
        params:
            ref= config['all']['REF'],
            RGinfo=getRGinfo
        threads: config['all']['THREADS']
        conda:
            "envs/bwa-mem.yaml"
        log: path.join(PATH_LOG, "bwa","{sample}.log")
        shell:
            "bwa mem -M -t {threads} -R {params.RGinfo} {params.ref} {input.fqtrim1} {input.fqtrim2} 2> {log} | samtools view -b -@ {threads}  -> {output} 2>> {log}"

        #"bwa aln -n {params.n} -t {threads} -q {params.q} {params.ref} {input} > {output.sai} 2> {log};"
        #"bwa samse -r '@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}' -f {output.samse}"
        #" {paramAbra=$JobName"_abra-realigned.bam"s.ref} {output.sai} {input} 2>> {log}"



rule Sambamba_sort:
    input:
        path.join(PATH_BAM, "{sample}_aligned_reads.bam"),
    output:
        qsorted=temp(path.join(PATH_BAM,"{sample}_querysorted.bam")),
        #sortbai=temp("../bam/{sample}_querysorted.bam.bai"),
    params:
        tmpdir=PATH_TEMP
    conda:
        "envs/sambamba.yaml"
    threads: config['all']['THREADS']
    log: path.join(PATH_LOG, "sambamba","{sample}_sort.log")
    shell:
       "sambamba sort -o {output.qsorted} --sort-by-name --tmpdir={params.tmpdir} -t {threads} -m 5000000000 {input} 2> {log} "

rule mark_duplicates:
    input:
        qsorted=path.join(PATH_BAM, "{sample}_querysorted.bam")
    output:
        dedup=temp(path.join(PATH_BAM, "{sample}_dedupped.bam")),
        metrics_file=path.join(PATH_QC, "Duplicates/{sample}.metrics.txt"),
    log: path.join(PATH_LOG, "mark_duplicates/{sample}.log"),
    params:
        tmpdir=PATH_TEMP
    threads: config['all']['THREADS'],
    conda:
        "envs/picard.yaml"
    shell:
        "picard -Xms4g -Xmx16g MarkDuplicates I={input.qsorted} O={output.dedup} M={output.metrics_file} "
        "ASSUME_SORT_ORDER=queryname TMP_DIR={params.tmpdir} &> {log} ;"

rule Sambamba_sort2:
    input:
        dedup=path.join(PATH_BAM,"{sample}_dedupped.bam"),
    output:
        coordsorted=path.join(PATH_BAM, "{sample}_coordsorted.bam"),
        sortbai=path.join(PATH_BAM, "{sample}_coordsorted.bam.bai"),
    params:
        tmpdir=PATH_TEMP,
    threads: config['all']['THREADS'],
    conda:
        "envs/picard.yaml"
    log: path.join(PATH_LOG, "sambamba/{sample}_sort2.log"),
    shell:
        "sambamba sort -o {output.coordsorted} --tmpdir {params.tmpdir} -t {threads} -m 5000000000 {input.dedup} 2> {log} && "
        "sambamba index -t {threads} {output.coordsorted} {output.sortbai} 2>> {log}"
        #"picard BuildBamIndex I={output.coordsorted} 2>> {log}"                      #change? samtools


rule CollectHsMetrics:
    input:
        coordsorted=path.join(PATH_BAM, "{sample}_coordsorted.bam"),
    output:
        HSmetrics=path.join(PATH_QC, "HSMetrics/{sample}_HSmetrics.txt"),
        PerTargetCov=path.join(PATH_QC, "HSMetrics/{sample}_PerTargetCov.txt"),
        PerBaseCov=path.join(PATH_QC, "HSMetrics/{sample}_PerBaseCov.txt")
    params:
        BAIT_INTERVALS=config["all"]["bait_intervals"],
        TARGET_INTERVALS=config["all"]["target_intervals"],
        ref=config["all"]["REF"],
    conda:
        "envs/picard.yaml"
    threads: config['all']['THREADS'],
    log: path.join(PATH_LOG, "picard/{sample}.metrics.log"),
    shell:
        "picard -Xmx16g CollectHsMetrics "
        "BAIT_INTERVALS={params.BAIT_INTERVALS} "
        "TARGET_INTERVALS={params.TARGET_INTERVALS} "
        "INPUT={input.coordsorted} "
        "OUTPUT={output.HSmetrics} "
        "METRIC_ACCUMULATION_LEVEL=ALL_READS "
        "PER_TARGET_COVERAGE={output.PerTargetCov} "
        "PER_BASE_COVERAGE={output.PerBaseCov} "
        "REFERENCE_SEQUENCE={params.ref} 2> {log}"


rule AlignmentMetrics:
    input:
        coordsorted=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        AlignmentMetrics=path.join(PATH_QC, "AlignmentMetrics/{sample}_AlignmentMetrics.txt")
    params:
        ref=config["all"]["REF"]
    conda:
        "envs/picard.yaml"
    log:
        path.join(PATH_LOG, 'picard/{sample}.alignmetrics.log')
    shell:
        "picard -Xmx16g CollectAlignmentSummaryMetrics "
        "R={params.ref} "
        "I={input.coordsorted} "
        "O={output.AlignmentMetrics} 2> {log}"

rule InsertsizeMetrics:
    input:
        coordsorted=path.join(PATH_BAM, "{sample}_coordsorted.bam"),
    output:
        InsertsizeMetrics=path.join(PATH_QC, "InsertsizeMetrics/{sample}_InsertsizeMetrics.txt"),
        Histogram=path.join(PATH_QC, "InsertsizeMetrics/{sample}_insertsize_histogram.pdf")
    conda:
        "envs/picard.yaml"
    log:
        path.join(PATH_LOG, 'picard/{sample}.insertsizemetrics.log')
    shell:
        "picard -Xmx16g CollectInsertSizeMetrics "
        "I={input.coordsorted} "
        "H={output.Histogram} "
        "O={output.InsertsizeMetrics} 2> {log}"


rule Metrics:
    input:
        HSmetrics=path.join(PATH_QC, "HSMetrics/{sample}_HSmetrics.txt"),
        AlignmentMetrics=path.join(PATH_QC, "AlignmentMetrics/{sample}_AlignmentMetrics.txt"),
        InsertsizeMetrics=path.join(PATH_QC, "InsertsizeMetrics/{sample}_InsertsizeMetrics.txt"),
    output:
        combined=path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt")
    script:
        path.join(PATH_PIPELINE, 'WES_Combine_metrics.R')

rule Collect_Metrics:
    input:
        expand(path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt"), sample = AllFiles.keys())
    output:
        path.join(PATH_QC, "Combined_Metrics.txt")
    run:
        out = list() 
        for f in input:
            out.append(pd.read_csv(f, sep=' '))
        out = pd.concat(out)
        out.to_csv(output[0], index=False)
            
#rule targetRegionFilter:
#    input:
#        coordsorted=path.join(PATH_BAM, "{sample}_coordsorted.bam")
#    output:
#        targetExonsBam=path.join(PATH_BAM, "{sample}_target_exons.bam"),
#    params:
#        targetExonsBed=config["all"]["target_bed"],
#    conda:
#        "envs/picard.yaml"
#        "envs/samtools.yaml"
#    log:
#        path.join(PATH_LOG, 'target_region_filter/{sample}.log')
#    shell:
#        "samtools view -L {params.targetExonsBed} {input.coordsorted} -b -o {output.targetExonsBam} 2> {log}"
#        "picard BuildBamIndex I={output.targetExonsBam} 2> {log}"
#        "java -Xmx8G -jar picard BuildBamIndex I={output.CapTargetBam} 2> {log}"

