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



##### Configurations VariantDetection ######
PATH_VAR = config['path']['variant']



## obtain sample list ###########
all_Samples=pd.read_csv(config['path']['sampleList'], sep='\t')

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
        path.join(PATH_QC, 'Combined_Metrics.txt'),
        #expand(path.join(PATH_QC, 'multiqc{trim}.html'), trim = ['', '_trim', '_bam']),
        #expand(path.join(PATH_QC, "ReadStats", "{sample}.reads.all"),
        #    sample = Tumor+Normal),
       #expand(path.join(PATH_VAR, 'LoFreq', '{sample}_somatic_final.snvs.vcf.gz'), sample=Tumor),
       # expand(path.join(PATH_VAR, "Mutect2/vcf/filtered/{sample}_Mutect2_passed.vcf.gz"), sample=Tumor)
       #expand("../vcf/filtered/{tumor}.all_evidenced.vcf", tumor=Tumor_samples.keys())
        path.join(PATH_VAR, 'funcotated/merged.vcf')

def getRGinfo(wildcards):
    fq=path.join(PATH_FASTQ, "trimmed", wildcards.sample+"_R1_trim.fq.gz")
    flowcell=subprocess.check_output("zcat "+fq+" | head -n 1  | awk 'BEGIN {FS=\":\"} {print $3}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()

    lane=subprocess.check_output("zcat "+fq+" | head -n 1 | awk 'BEGIN { FS=\":\"} {print $4}' ", preexec_fn=lambda:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL), shell=True).decode('utf-8').strip()
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
    RGinfo="\'@RG\tID:"+lane+"\tSM:"+name+"\tPL:ILLUMINA\tLB:SeqCap\tPU:"+flowcell+"\'"

    return(RGinfo)

def getNormalSample(wildcards):
    normal = pairs[wildcards.sample]
    name=re.match('[a-zA-Z0-9\-]*', normal).group(0)
    return(name)

def getTumorSample(wildcards):
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
    return(name)



if PLATFORM in ['PE', 'pe']: # I only implement PE for now

    rule merge_fastq:
        input:
            lambda wildcards: AllFiles[wildcards.sample][wildcards.R1R2]
        output:
            path.join(PATH_FASTQ, '{sample}.fastq.merged{R1R2}.gz')
        threads: 10
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
        threads: 1 #config['all']['THREADS']
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
            html1=path.join(PATH_QC,"fastqc", "trimmed","{sample}_R1_trim_fastqc.html"),
            html2=path.join(PATH_QC,"fastqc", "trimmed","{sample}_R2_trim_fastqc.html"),
            zip1=temp(path.join(PATH_QC, "fastqc", "trimmed", "{sample}_R1_trim_fastqc.zip")),
            zip2=temp(path.join(PATH_QC, "fastqc", "trimmed", "{sample}_R2_trim_fastqc.zip")),
        threads: 1 #config['all']['THREADS']
        conda:
            "envs/fastqc.yaml"
        params:
            fastqc_dir=path.join(PATH_QC, "fastqc", "trimmed")
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
            expand(path.join(PATH_QC, "fastqc", "trimmed", "{sample}{pair}_trim_fastqc.zip"),
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


rule ReadsStats:
    input:
        bam=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        path.join(PATH_QC, "ReadStats", "{sample}.reads.all")
    params:
        outdir = path.join(PATH_QC, "ReadStats"),
        path_snakemake = srcdir('.')
    conda:
            "envs/bwa-mem.yaml"
    shell:
        """
        {params.path_snakemake}/scripts/stats.sh {input.bam} {wildcards.sample} {params.outdir} {output}
        """

rule fastqc_bam:
    input:
        bam=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        bamqc = path.join(PATH_QC, "fastqc", "bam", "{sample}_coordsorted_fastqc.html"),
        bamqczip = path.join(PATH_QC, "fastqc", "bam", "{sample}_coordsorted_fastqc.zip")
    threads: 1 #config['all']['THREADS']
    params:
        fastqc_dir=path.join(PATH_QC, "fastqc", "bam")
    conda:
        "envs/fastqc.yaml"
    log: path.join(PATH_LOG, "fastqc","{sample}.bam.log")
    shell:
        "fastqc {input.bam} --bam_mapped --outdir {params.fastqc_dir} --threads {threads} 2> {log}"  #Java??



rule multiqc_bam:
        input:
            expand(path.join(PATH_QC, "fastqc", "bam", "{sample}_coordsorted_fastqc.zip"),
                            sample=AllFiles.keys())
        output:
            path.join(PATH_QC,'multiqc_bam.html') # empty wildcard allowed
        log:
            path.join(PATH_LOG,'multiqc_bam.log')
        wrapper:
            '0.35.0/bio/multiqc'



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
    threads:  1 #config['all']['THREADS'],
    log: path.join(PATH_LOG, "picard/{sample}.metrics.log"),
    shell:
        "picard -Xmx16g CollectHsMetrics "
        "BAIT_INTERVALS={params.BAIT_INTERVALS} "
        "TARGET_INTERVALS={params.TARGET_INTERVALS} "
        "COVERAGE_CAP=10000 " # to go beyond 200 reads
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
        AlignmentMetrics=path.join(PATH_QC, "AlignmentMetrics","{sample}_AlignmentMetrics.txt")
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
        """picard -Xmx16g CollectInsertSizeMetrics I={input.coordsorted} H={output.Histogram} O={output.InsertsizeMetrics} 2> {log} 
	touch -a {output.InsertsizeMetrics}
	touch -a {output.Histogram}
	"""

#        [[ ! -f "{output.InsertsizeMetrics}" ]] && echo "" > {output.InsertsizeMetrics} fi
#        [[ ! -f "{output.Histogram}" ]] && echo "" > {output.Histogram} fi"""


rule Metrics:
    input:
        HSmetrics=path.join(PATH_QC, "HSMetrics/{sample}_HSmetrics.txt"),
        AlignmentMetrics=path.join(PATH_QC, "AlignmentMetrics/{sample}_AlignmentMetrics.txt"),
        InsertsizeMetrics=path.join(PATH_QC, "InsertsizeMetrics/{sample}_InsertsizeMetrics.txt"),
        readStat=path.join(PATH_QC, "ReadStats", "{sample}.reads.all")
    output:
        combined=path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt")
    script:
        path.join(PATH_PIPELINE, 'scripts', 'WES_Combine_metrics.R')

rule Collect_Metrics:
    input:
        expand(path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt"), sample = AllFiles.keys())
    output:
        path.join(PATH_QC, "Combined_Metrics.txt")
    params:
        columns = ['SAMPLE', 'TOTALREADS', 'Q37_READS', 'PF_UNIQUE_READS', 'PCT_PF_UQ_READS', 'MEAN_TARGET_COVERAGE',
                'MEDIAN_TARGET_COVERAGE', 'PCT_USABLE_BASES_ON_TARGET', 'PCT_TARGET_BASES_30X', 'PCT_TARGET_BASES_100X', 'PCT_CHIMERAS',
                'MEAN_READ_LENGTH', 'MEDIAN_INSERT_SIZE', 'PCT_EXC_DUPE', 'PCT_EXC_OVERLAP', 'PCT_EXC_OFF_TARGET']
    run:
        out = list() 
        for f in input:
            out.append(pd.read_csv(f, sep=' '))
        out = pd.concat(out)
        out = out[params.columns]
        out.to_csv(output[0], index=False)
            


############## VARIANT DETECTION #################


rule LoFreq:
    input:
        normal = lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample] + '_coordsorted.bam'),
        tumor=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        snps=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.snvs.vcf.gz"),
        indels=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.indels.vcf.gz")
    params:
        ref=config["all"]["REF"],
        threads=config["all"]["THREADS"],
        min_cov=config["LoFreq"]["min_cov"],
        min_mq=config["LoFreq"]["min_mq"],
        min_bq=config["LoFreq"]["min_bq"],
        min_alt_bq=config["LoFreq"]["min_alt_bq"],
        max_depth=config["LoFreq"]["max_depth"],
        sig=config["LoFreq"]["sig"],
        prefix=path.join(PATH_VAR, "LoFreq", "{sample}_")
    conda:
        path.join(PATH_PIPELINE, "envs", "lofreq.yaml")
    log: path.join(PATH_LOG, "LoFreq", "{sample}_lofreq.txt")
    shell:
        """
        lofreq somatic \
          -n {input.normal} -t {input.tumor} -f {params.ref} -o {params.prefix} \
        --min-cov {params.min_cov} \
        --threads {params.threads} --call-indels --verbose \
        2> {log}
        """

rule LoFreq_combine:
    input:
        snps=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.snvs.vcf.gz"),
        indels=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.indels.vcf.gz")
    output:
        tmp=temp(path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.combined.tmp.vcf")),
        out=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.combined.vcf")
    conda:
        "envs/lofreq.yaml"
    shell:
        """
	tabix -f {input.snps} &&
	tabix -f {input.indels} &&
        vcfcombine {input.snps} {input.indels} > {output.tmp} && 
        sed 's/\%//' {output.tmp} | sed 's/FREQ\,Number\=1\,Type\=String/FREQ\,Number\=1\,Type\=Float/' > {output.out}
        """


rule LoFreq_readStatFilter:
    input:
        raw_vcf=path.join(PATH_VAR, "LoFreq", "{sample}_somatic_final.combined.vcf")
    output:
        tmp_vcf=temp(path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_tmp.vcf")),
        unsorted_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_unsorted.vcf"),
        filtered_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_filt.vcf"),
    params:
        af_min=config["LoFreq"]["af_min"],
        cov_min=config["LoFreq"]["cov_min"],
        sb_alpha=config["LoFreq"]["sb_alpha"],
        #SnpSift_filter=config["LoFreq"]["SnpSift_filter"],
        hg19_dict=config["all"]["HG19_DICT"],
    conda:
        "envs/lofreq.yaml"
    log: path.join(PATH_LOG, "LoFreq/{sample}_lofreq_readStatFilter.txt")
    shell:
        """
        lofreq filter --verbose --af-min {params.af_min} --cov-min {params.cov_min} --sb-alpha {params.sb_alpha} --sb-incl-indels -i {input.raw_vcf} -o {output.tmp_vcf}
        SnpSift filter -f {output.tmp_vcf} "(DP4[2]>2) & (DP4[3]>2) & ((na HRUN) | (HRUN<8))" > {output.unsorted_vcf}
        picard SortVcf I={output.unsorted_vcf} O={output.filtered_vcf} SD={params.hg19_dict} 2> {log}
        """

rule zip_vcf:
    input:
        path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf')
    output:
        vcf=path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf.gz'),
        idx=path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf.gz.tbi')
    conda:
        "envs/lofreq.yaml"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule Mutect2:
    input:
        normal=lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample] + '_coordsorted.bam'),
        tumor=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        raw=path.join(PATH_VAR, "Mutect2/{sample}_Mutect2.vcf"),
        f1r2=path.join(PATH_VAR, "Mutect2/{sample}_f1r2.tar.gz")
    params:
        ref=config["all"]["REF"],
        pon=config["Mutect2"]["pon"],
        gnomad=config["Mutect2"]["gnomad"],
        interval=config["Mutect2"]["interval"],
        normal=getNormalSample,
        tumor=getTumorSample
    conda: "envs/gatk4.yaml"
    threads: config["all"]["THREADS"]
    log:	path.join(PATH_LOG, "Mutect2/{sample}_mutect2.txt")
    shell:
        """
        gatk Mutect2 -R {params.ref} -I {input.normal} -I {input.tumor} \
        -O {output.raw} --normal-sample {params.normal} --tumor-sample {params.tumor} \
        --germline-resource {params.gnomad} \
        -pon {params.pon} -L {params.interval} --f1r2-tar-gz {output.f1r2} --max-mnp-distance 0 --native-pair-hmm-threads {threads} 2> {log}
        """


#pass this raw mutect data to LearnReadOrientationModel - to be able to filter on orientation bias
rule LearnReadOrientationModel:
    input:
        f1r2=path.join(PATH_VAR, "Mutect2/{sample}_f1r2.tar.gz")
    output:
        obmodel=path.join(PATH_VAR, "Mutect2/{sample}_read-orientation-model.tar.gz")
    conda:
        "envs/gatk4.yaml"
    log:
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2_orientation.txt")
    shell:
        """
        gatk LearnReadOrientationModel -I {input.f1r2} -O {output.obmodel} 2> {log}
        """

# this is to determine contamination
rule CalculateContamination:
    input:
        exonbam=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        pileup=temp(path.join(PATH_VAR, "Mutect2/{sample}_getpileupsummaries.table")),
        segments=path.join(PATH_VAR, "Mutect2/{sample}_segments.table"),
        contamination=path.join(PATH_VAR, "Mutect2/{sample}_calculatecontamination.table")
    params:
        variants=config["Mutect2"]["variants"],
        interval=config["Mutect2"]["interval"],
    conda:
        "envs/gatk4.yaml"
    log: 
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2.txt")
    shell:
        "gatk GetPileupSummaries -I {input.exonbam} -O {output.pileup} -V {params.variants} -L {params.interval} 2>> {log} &&"
        "gatk CalculateContamination -I {output.pileup} -tumor-segmentation {output.segments} -O {output.contamination} 2>> {log}"


#pass the learned read orientation model and contaminationtable in FilterMutectCalls
rule FilterMutectCalls:
    input: 
        raw=path.join(PATH_VAR, "Mutect2/{sample}_Mutect2.vcf"),
        obmodel=path.join(PATH_VAR, "Mutect2/{sample}_read-orientation-model.tar.gz"),
        segments=path.join(PATH_VAR, "Mutect2/{sample}_segments.table"),
        contamination=path.join(PATH_VAR, "Mutect2/{sample}_calculatecontamination.table"),
    output:
        filtered=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_filt.vcf"),
    params:
        ref=config["all"]["REF"],
        events=config["Mutect2"]["clustered"],      #50		#this is for the clustered events, default =2 but than artifacts cause that real events are filtered and SHM is filtered out
        af_min=config["Mutect2"]["af_min"],
        reads= config["Mutect2"]["reads_per_strand"],	
        min_reads= config["Mutect2"]["min_reads"],	
    conda:
        "envs/gatk4.yaml"
    log: 
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2.txt")
    shell:
        """
        gatk FilterMutectCalls -V {input.raw} -R {params.ref} \
        --max-events-in-region {params.events} --min-allele-fraction {params.af_min} \
        --tumor-segmentation {input.segments} --contamination-table {input.contamination} \
        --min-reads-per-strand {params.reads} --unique-alt-read-count {params.min_reads} \
        --ob-priors {input.obmodel} -O {output.filtered} 2>> {log} 
        """

rule Mutect_passed: 
    input:
        filtered=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_filt.vcf"),
    output:
        passed=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.vcf"),
        pass_csv=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.csv"),
        passed_gz=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.vcf.gz"),
    params:
        fields='CHROM POS REF ALT DP "GEN[0].AF" "GEN[0].F1R2" "GEN[0].F2R1" FILTER',
    conda: 
        "envs/SnpSift.yaml"
    shell:
        """
        SnpSift filter -f {input.filtered} "(FILTER = 'PASS' & (ROQ > 20 | TLOD > 20))" > {output.passed} &&
        SnpSift extractFields -e "." {output.passed} {params.fields}> {output.pass_csv} &&
        pbgzip -c {output.passed} > {output.passed_gz} &&
        tabix -s1 -b2 -e2 {output.passed_gz}
        """


rule Intersect_VariantCalls:
    input:
        lofreq_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_filt.vcf.gz"),
        Mutect2_vcf=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.vcf.gz"),
    output:
        intersect_vcf=path.join(PATH_VAR, "intersect/{sample}_intersect.vcf"),
        outersect_vcf=path.join(PATH_VAR, "intersect/{sample}_outersect.vcf")
    conda:
        "envs/SnpSift.yaml"
    log: path.join(PATH_LOG, "Intersect_variantCalls/{sample}_intersectCalls.txt")  
    shell:
        """
        bcftools isec -n=2 -w1 {input.Mutect2_vcf} {input.lofreq_vcf} -o {output.intersect_vcf} -O v 2> {log}
        bcftools isec {input.Mutect2_vcf} {input.lofreq_vcf} -c all -n +0 -o {output.outersect_vcf} -O v 2>> {log}
        """

rule Annotate_VariantCalls:
    input:
        intersect_vcf=path.join(PATH_VAR, "intersect/{sample}_intersect.vcf")
    output:
        annotated=path.join(PATH_VAR, "annotated/{sample}.annotated.vcf")
    params:
        Java_mem=config["all"]["Java_mem"],
        dbsnp=config["all"]["dbsnp"],
        clinvar=config["all"]["clinvar"],
        Cosmic=config["all"]["Cosmic"],
        gnomAD=config["all"]["gnomAD"],
        HMF_PON=config["all"]["HMF_PON"],
    conda:
        "envs/SnpSift.yaml"
    log:
        path.join(PATH_LOG, "Annotate_variantCalls/{sample}_annotation.txt")
    shell:
        """
        SnpSift annotate {params.Java_mem} -v {params.Cosmic} {input.intersect_vcf} 2>> {log} |
        SnpSift annotate {params.Java_mem} -v -info 'gnomAD_AF' {params.gnomAD} - 2>> {log} |
        SnpSift annotate {params.Java_mem} -v -info 'PON_COUNT' {params.HMF_PON} - > {output.annotated} 2>> {log}
        """


rule Funcotator:
    input:
        annotated=path.join(PATH_VAR, "annotated/{sample}.annotated.vcf"),
    output:
        funcotated=path.join(PATH_VAR, "funcotated/{sample}_funcotated.vcf"),
    params:
        datasource=config["Funcotator"]["datasource"],
        ref=config["all"]["REF"],
        version=config["Funcotator"]["hg_version"],
        file_format=config["Funcotator"]["file_format"],
#        transcriptlist=config["Funcotator"]["transcriptlist"],
    conda:
        "envs/gatk4.yaml"
    log:	
        path.join(PATH_LOG, "Funcotator/{sample}_funcotator.txt")
    shell:
        """
        gatk Funcotator \
     		--variant {input.annotated}\
        	--reference {params.ref} \
     		--ref-version {params.version} \
        	--data-sources-path {params.datasource} \
     		--output {output.funcotated} \
     		--output-file-format {params.file_format} 2>> {log}
        """

rule merge_vcf:
    input:
        expand(path.join(PATH_VAR, 'funcotated/{sample}_funcotated.vcf'), sample=Tumor)
    output:
        path.join(PATH_VAR, 'funcotated/merged.vcf')
    conda:
        "envs/lofreq.yaml"
    log: path.join(PATH_LOG, "merge_vcfs/funcotated.txt")  
    shell:
        """
        vcfcombine {input} > {output} 2> {log}
        """



rule ExtractVcfFields_csv:
    input:
        funcotated=path.join(PATH_VAR, "funcotated/{sample}_funcotated.vcf"),
    output:
        tmp_csv_snpsift=temp(path.join(PATH_VAR, "funcotated/{sample}_ann.csv")),
        tmp_csv_func=temp(path.join(PATH_VAR, "funcotated/{sample}_func.csv")),
        funcotated_csv=temp(path.join(PATH_VAR, "funcotated/{sample}_funcotated_org_head.csv")),
    params:
        fields='CHROM POS REF ALT DP "GEN[0].F1R2" "GEN[0].F2R1" "GEN[0].AF" \
        AS_SB_TABLE "GEN[0].SB" GERMQ MBQ MFRL MMQ MPOS ROQ RPA RU STRQ STR SEQQ \
            STRANDQ  TLOD COSM.ID AA CDS FATHMM MUT.ST SNP PON PON_COUNT POPAF gnomAD_AF',
    conda:
        "envs/SnpSift.yaml"
    shell:
        """
        SnpSift extractFields -e "." {input.funcotated} {params.fields} > {output.tmp_csv_snpsift} &&
        ./scripts/createFuncotationCsvFromFuncotatorVcf.sh {input.funcotated} > {output.tmp_csv_func} &&
        paste {output.tmp_csv_snpsift} {output.tmp_csv_func} > {output.funcotated_csv}
        """		



	#TODO when a table is empty an error occurs and no output is created - maybe create if? 
	#TODO PON_COUNT < 4 is not recognized correctly, it seems 25 is counted as 2.5 or something like that. if -I is not used F1R2 is not recognized as comma separated list but as decimal value. quick fix is to filter PON_COUNT > 4 out in mutation_overviews.
# (SNP is null OR SNP = 'false') AND - removed this, because its from the Cosmic db and less reliable


