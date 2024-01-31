#/bin/bash
#import os, sys
import subprocess
import pandas as pd
import glob
import os.path as path
configfile: "config.yaml"
#import os
#import sys
#os.system("check_config")
#import shlex
#output = subprocess.call(['check_config'])
report: "report/workflow.rst"

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
#contains all the tumor nonpaired samples in a list
new_non_paired = []
#contain all the tumor id of the paired samples
new_Tumor_Matched = []
# contains all the pairs of tumor-matched + normal. where the Tumor is the key. 
new_pairs = {}
#contains all tumor sample id
new_all_list = []

for tumor,normal in zip(all_Samples['Tumor'], all_Samples['Normal']):
    print(tumor,normal)
    if isinstance(normal, float): #means normal = nan
        new_non_paired += [tumor]
        new_all_list += [tumor]
    else:
        new_pairs[tumor] = normal
        new_Tumor_Matched += [tumor]
        new_all_list += [tumor]
        #new_all_list += [normal]

#If statements can be build in if you have sample list with multiple dublicates. 
#        if normal not in all_samples:
#           all_samples += [normal]
#    if tumor not in all_samples:
#           all_samples += [tumor]

## obtain paired/tumorOnly sample lists ###########
paired_Samples= all_Samples[~all_Samples['Normal'].isnull()]
tumorOnly_Samples =  all_Samples[all_Samples['Normal'].isnull()]

Tumor=list(all_Samples['Tumor']) # these are IDs, and use as the key to fetch files.
Normal=list(all_Samples['Normal'])


def getnames(samplelist, platform, pathdata):
    Files = []
    RNAIDs = []
    SAMPLES = dict()
    for sample in samplelist:
        if str(sample) == 'nan':
            pass
        else:
            if platform in ['SR', 'sr']:
                for prefix in PREFIX:
                    SAMPLES[sample] = glob.glob(path.join(pathdata, sample+'*'+prefix+'.fastq.gz'))
                    SAMPLES[sample].sort()
            if platform in ['PE', 'pe']:
                SAMPLES[sample] = dict()
                SAMPLES[sample]['R1'] = glob.glob(path.join(pathdata, str(sample)+'*'+ PREFIX[0] +'.fastq.gz'))   #if sample is NA thinks it is a float, fixed by str()
                SAMPLES[sample]['R1'].sort()
                SAMPLES[sample]['R2'] = glob.glob(path.join(pathdata, str(sample)+'*'+PREFIX[1]+'.fastq.gz'))       #if sample is NA thinks it is a float, fixed by str()
                SAMPLES[sample]['R2'].sort()
    return(SAMPLES)

Tumor_samples = getnames(Tumor, PLATFORM, PATH_FASTQ)
Normal_samples = getnames(Normal, PLATFORM, PATH_FASTQ)
AllFiles =  getnames(Tumor, PLATFORM, PATH_FASTQ) #Tumor_samples  changed this, gave problem with tumor sample. 
AllFiles.update(Normal_samples)
Paired_samples_Id = list(Normal_samples) + (list(Tumor_samples.keys())[0:len(Normal_samples)])


Non_Paired = list(Tumor_samples.keys())[len(Normal_samples):]
Tumor_Matched = (list(Tumor_samples.keys())[0:len(Normal_samples)])
pairs = dict(zip(paired_Samples['Tumor'],paired_Samples['Normal']))
all_list = list(Non_Paired) + list(Tumor_Matched)

#contains the location of all fastq.gz file for all samples.
#Is an dictionary with key the sample id, and as value another dictionary with key R1 and R2 that contain a list with the location of the R1/R2 fastq.gz files. 
print('ALLFILES', AllFiles)


#New build
Non_Paired = new_non_paired
Tumor_Matched = new_Tumor_Matched
pairs = new_pairs
all_list = new_all_list


#####################################


rule all:
    input:
        expand(path.join(PATH_VAR, "Mutect2_paired/{sample_paired}_Mutect2.vcf"), sample_paired = Tumor_Matched), 			#paired_Samples['Tumor']),
        expand(path.join(PATH_VAR, "Mutect2_tumorOnly/{sample_tumor}_Mutect2.vcf"), sample_tumor = Non_Paired), 			#tumorOnly_Samples['Tumor']),        
        expand(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_final.snvs.vcf.gz"), sample_paired = Tumor_Matched), 	#paired_Samples['Tumor']),  
        expand(path.join(PATH_VAR, "LoFreq_tumorOnly", "{sample_tumor}_somatic_final.combined.vcf"), sample_tumor = Non_Paired), 	#tumorOnly_Samples['Tumor']),         
        expand(path.join(PATH_VAR, "Mutect2/{sample}_read-orientation-model.tar.gz"), sample = all_list),  
        expand(path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_filt.vcf"),  sample = all_list),
        path.join(PATH_QC, 'Combined_Metrics.txt'),
        expand(path.join(PATH_QC, 'multiqc{trim}.html'), trim = ['', '_trim', '_bam']),  
        path.join(PATH_VAR, 'funcotated/merged.vcf'),
        path.join(PATH_VAR, 'merged_mutation_counts.tsv'), 
        expand(path.join(PATH_VAR, "mafs/{sample}.maf"),sample = all_list),   			 #Might not be needed anymore
        path.join(PATH_VAR, 'oncoplot.pdf'),							 #new added
        expand("benchmarks/merged/{sample}_merged_benchmark.tsv", sample=AllFiles.keys()),	 #new added
       


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
    normal = pairs[wildcards.sample_paired]
    name=re.match('[a-zA-Z0-9\-]*', normal).group(0)
    return(name)

def getTumorSample(wildcards):
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample_paired).group(0)
    return(name)

if PLATFORM in ['PE', 'pe']: #only PE is implemented

    #THIS RULE IS ONLY NEEDED IF NO FASTQ FILES ARE PRESENT
    #rule bam_to_fastq:
    #    input: 
    #        path.join("/net/beegfs/cfg/tgac/ferdinand/barabara_wes_25_09_2023/fake", "{sample}.bam")  
    #    output:
    #        fq1=path.join(PATH_FASTQ,  "{sample}_R1_001.fastq.gz"),
    #        fq2=path.join(PATH_FASTQ, "{sample}_R2_001.fastq.gz"), 
    #        sort=temp(path.join(PATH_FASTQ, "{sample}_sorted.bam")),
    #    conda:
    #        "envs/samtools2.yaml"
    #    shell:
    #        """
    #        samtools sort -n {input} -o {output.sort} 
    #        samtools fastq -@ 8 {output.sort} -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n 
    #        """

    rule merge_fastq: 
        input:
            input= (lambda wildcards: AllFiles[wildcards.sample][wildcards.R1R2]),
            fake=  (lambda wildcards:  expand(path.join(PATH_FASTQ, "{sample}_sorted.bam"), sample=AllFiles.keys())),
        output:
            path.join(PATH_FASTQ, '{sample}.fastq.merged{R1R2}.gz')
        #threads: 10  #might not need treads
        shell:
            """
            test -f {input}  && echo True || echo Error  
            cat {input.input} > {output}
            """
    rule fastqc:
        input:
            fq1=path.join(PATH_FASTQ, "{sample}.fastq.mergedR1.gz"),
            fq2=path.join(PATH_FASTQ, "{sample}.fastq.mergedR2.gz")
        output:
            html1=path.join(PATH_QC, "fastqc", "raw", "{sample}.fastq.mergedR1_fastqc.html"),
            html2=path.join(PATH_QC, "fastqc", "raw", "{sample}.fastq.mergedR2_fastqc.html"),
            zip1=temp(path.join(PATH_QC, "fastqc", "raw", "{sample}.fastq.mergedR1_fastqc.zip")),
            zip2=temp(path.join(PATH_QC, "fastqc", "raw", "{sample}.fastq.mergedR2_fastqc.zip"))
        params:
            fastqc_dir=path.join(PATH_QC, "fastqc", "raw")
        conda:
            "envs/fastqc.yaml"
        threads:  config['all']['THREADS'],
        benchmark:
            "benchmarks/{sample}/rule_fastqc.bwa.benchmark.txt"
        log: path.join(PATH_LOG, "fastqc","{sample}.log")
        shell:
            "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"  

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
            adapter5=config['cutadapt']['adapter'][1],   #standard ILLUMINA but also IDT
            min_len=config["trim"]["min_len"],
        log: path.join(PATH_LOG, "trimming", "{sample}.log"),
        benchmark:
            "benchmarks/{sample}/rule_cutadapt_trimming.bwa.benchmark.txt"
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
        conda:
            "envs/fastqc.yaml"
        params:
            fastqc_dir=path.join(PATH_QC, "fastqc", "trimmed")
        log: path.join(PATH_LOG, "fastqc","{sample}.trim.log")
        benchmark:
            "benchmarks/{sample}/fastqc_trimmed.bwa.benchmark.txt"
        shell:
            "fastqc {input.fq1} {input.fq2} --outdir {params.fastqc_dir} --threads {threads} 2> {log}"

   
    rule new_multiqc_raw:
        input:
            expand(path.join(PATH_QC, "fastqc", "raw", "{sample}{pair}_fastqc.zip"), pair=['.fastq.mergedR1','.fastq.mergedR2'], sample=AllFiles.keys()) 
        output:
            report(path.join(PATH_QC, 'multiqc.html'), category="Multiqc")  
        log:
            path.join(PATH_LOG,'multiqc.log')
        params:
            output_dir=path.join(PATH_QC)
        conda: "envs/multiqc_env_new.yaml"
        shell:
            """
            multiqc {input} -o {params.output_dir}  -n "multiqc.html" --force
            """
    #wrapper problem new rule made called "new_multiqc_raw"       error was "dataclasses not found"
    #rule multiqc_raw:
    #    input:
    #        expand(path.join(PATH_QC, "fastqc", "raw", "{sample}{pair}_fastqc.zip"), pair=['.fastq.mergedR1','.fastq.mergedR2'], sample=AllFiles.keys()) #added this line. 
    #    output:
    #        report(path.join(PATH_QC, 'multiqc.html'), category="Multiqc")   #empty wildcard allowed
    #    log:
    #        path.join(PATH_LOG,'multiqc.log')
    #    wrapper:
    #         "v1.25.0/bio/multiqc"  #'0.35.0/bio/multiqc'   
   

    rule new_multiqc_trim:
        input:
            expand(path.join(PATH_QC, "fastqc", "trimmed", "{sample}{pair}_trim_fastqc.zip"), pair=['_R1','_R2'], sample=AllFiles.keys())
        output:
            report(path.join(PATH_QC,'multiqc_trim.html'), category="Multiqc") 
        log:
            path.join(PATH_LOG,'multiqc_trim.log')
        params:
            output_dir=path.join(PATH_QC)
        conda: "envs/multiqc_env_new.yaml"
        shell:
            """
            multiqc {input} -o {params.output_dir}  -n "multiqc_trim.html" --force
            """

    #wapper problem new rule made called "new_multiqc_trim"      error was: "dataclasses not found"
    #rule multiqc_trim:
    #    input:
    #       expand(path.join(PATH_QC, "fastqc", "trimmed", "{sample}{pair}_trim_fastqc.zip"), pair=['_R1','_R2'], sample=AllFiles.keys())
    #    output:
    #        report(path.join(PATH_QC,'multiqc_trim.html'), category="Multiqc") # empty wildcard allowed 
    #    log:
    #        path.join(PATH_LOG,'multiqc_trim.log')
    #    wrapper:
    #         "v1.25.0/bio/multiqc" #'0.35.0/bio/multiqc'       

    rule bwa_mem:
        input:
            fqtrim1=path.join(PATH_FASTQ, "trimmed", "{sample}_R1_trim.fq.gz"),
            fqtrim2=path.join(PATH_FASTQ, "trimmed", "{sample}_R2_trim.fq.gz"),
        output:
             temp(path.join(PATH_BAM, "{sample}_aligned_reads.bam")),
        params:
            ref= config['all']['REF_bwa'],
            RGinfo=getRGinfo                       #this is calling the function getRGinfo(wildcards)
        threads: config['all']['THREADS']
        resources:
            mem_mb = 17000
        conda:
            "envs/bwa-mem.yaml"
        log: path.join(PATH_LOG, "bwa","{sample}.log") 
        benchmark:
            "benchmarks/{sample}/bwa_mem.bwa.benchmark.txt"
        shell:
            "bwa mem -M -t {threads} -R {params.RGinfo} {params.ref} {input.fqtrim1} {input.fqtrim2} 2> {log} | samtools view -b -@ {threads}  -> {output} 2>> {log}"


rule Sambamba_sort:
    input:
        path.join(PATH_BAM, "{sample}_aligned_reads.bam"),
    output:
        qsorted=temp(path.join(PATH_BAM,"{sample}_querysorted.bam"))
    params:
        tmpdir=PATH_TEMP
    conda:
        "envs/sambamba.yaml"
    threads: config['all']['THREADS']
    log: path.join(PATH_LOG, "sambamba","{sample}_sort.log")
    benchmark:
            "benchmarks/{sample}/sambama_sort.bwa.benchmark.txt"
    shell:
        "sambamba sort -o {output.qsorted} --sort-by-name --tmpdir={params.tmpdir} -t {threads} -m 5000000000 {input} 2> {log} "

rule mark_duplicates:
    input:
        qsorted=path.join(PATH_BAM, "{sample}_querysorted.bam")
    output:
        dedup=temp(path.join(PATH_BAM, "{sample}_dedupped.bam")),
        metrics_file=path.join(PATH_QC, "Duplicates/{sample}.metrics.txt")
    log: path.join(PATH_LOG, "mark_duplicates/{sample}.log")
    params:
        tmpdir=PATH_TEMP
    threads: 1           #can be increased if input file is not to large.
    conda:
        "envs/picard.yaml"
    benchmark:
            "benchmarks/{sample}/mark_duplicates.bwa.benchmark.txt"
    resources:
        mem_mb = 400000, #can be decreased if input file is small.   		#lambda wildcards, input, attempt: (input.size//1000000) * attempt * 5,
        disk_mb= 400000, #can be decreased if input file is small.    		#lambda wildcards, input, attempt: (input.size//1000000) * attempt * 5, 
    shell:
       "picard -Xms360g -Xmx360g MarkDuplicates COMPRESSION_LEVEL=9 sorting_collection_size_ratio=0.0001  I={input.qsorted} O={output.dedup} M={output.metrics_file} "
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
    resources:
        mem_mb = lambda wildcards, input, attempt: (input.size//1000000) * attempt * 3, 	#can also be tested with 400000,
    conda:
        "envs/picard.yaml"
    log: path.join(PATH_LOG, "sambamba/{sample}_sort2.log"),
    benchmark:
            "benchmarks/{sample}/sambamba_sort2.bwa.benchmark.txt"
    shell:
        "sambamba sort -o {output.coordsorted} --tmpdir {params.tmpdir} -t {threads} -m 5000000000 {input.dedup} 2> {log} && "
        "sambamba index -t {threads} {output.coordsorted} {output.sortbai} 2>> {log}"

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
    benchmark:
            "benchmarks/{sample}/readsStats.bwa.benchmark.txt"
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
    params:
        fastqc_dir=path.join(PATH_QC, "fastqc", "bam")
    conda:
        "envs/fastqc.yaml"
    log: path.join(PATH_LOG, "fastqc","{sample}.bam.log")
    benchmark:
            "benchmarks/{sample}/fastqc_bam.bwa.benchmark.txt"
    shell:
        "fastqc {input.bam} --bam_mapped --outdir {params.fastqc_dir} --threads {threads} 2> {log}"  


rule new_multiqc_bam:
    input:
        expand(path.join(PATH_QC, "fastqc", "bam", "{sample}_coordsorted_fastqc.zip"),
                            sample=AllFiles.keys())
    output:
        report(path.join(PATH_QC,'multiqc_bam.html'), category="Multiqc") # empty wildcard allowed               
    log:
        path.join(PATH_LOG,'multiqc_bam.log')
    params:
            output_dir=path.join(PATH_QC)
    conda: "envs/multiqc_env_new.yaml"
    shell:
            """
            multiqc {input} -o {params.output_dir}  -n "multiqc_bam.html" --force
            """

#problems with wrapper new rule made called "new_multiqc_bam"
#rule multiqc_bam:
#    input:
#        expand(path.join(PATH_QC, "fastqc", "bam", "{sample}_coordsorted_fastqc.zip"),
#                            sample=AllFiles.keys())
#    output:
#        report(path.join(PATH_QC,'multiqc_bam.html'), category="Multiqc") # empty wildcard allowed                   
#    log:
#        path.join(PATH_LOG,'multiqc_bam.log') 
#    wrapper:
#        'v1.25.0/bio/multiqc' #'0.35.0/bio/multiqc'


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
    threads:  1 		#can also be changed to config['all']['THREADS'],
    resources:
        mem_mb = 17000
    log: path.join(PATH_LOG, "picard/{sample}.metrics.log"),
    benchmark:
            "benchmarks/{sample}/collectHsMetrics.bwa.benchmark.txt"
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
    resources:
        mem_mb = 17000
    log:
        path.join(PATH_LOG, 'picard/{sample}.alignmetrics.log')
    benchmark:
            "benchmarks/{sample}/alignmentMetrics.bwa.benchmark.txt"
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
        Histogram=report(path.join(PATH_QC, "InsertsizeMetrics/{sample}_insertsize_histogram.pdf"), category="Metrics")   
    conda:
        "envs/picard.yaml"
    log:
        path.join(PATH_LOG, 'picard/{sample}.insertsizemetrics.log')
    resources:
        mem_mb = 17000
    benchmark:
            "benchmarks/{sample}/insertsizeMetrics.bwa.benchmark.txt"
    shell:
        """picard -Xmx16g CollectInsertSizeMetrics I={input.coordsorted} H={output.Histogram} O={output.InsertsizeMetrics} 2> {log} 
	touch -a {output.InsertsizeMetrics}
	touch -a {output.Histogram}
	"""

rule Metrics:
    input:
        HSmetrics=path.join(PATH_QC, "HSMetrics/{sample}_HSmetrics.txt"),
        AlignmentMetrics=path.join(PATH_QC, "AlignmentMetrics/{sample}_AlignmentMetrics.txt"),
        InsertsizeMetrics=path.join(PATH_QC, "InsertsizeMetrics/{sample}_InsertsizeMetrics.txt"),
        readStat=path.join(PATH_QC, "ReadStats", "{sample}.reads.all"),
        #check=path.join(PATH_QC, "ReadStats", "{sample}.reads.q37")
        #check=path.join(PATH_QC, "ReadStats", "{sample}.reads.aligned") 
    output:
        combined=path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt")
    benchmark:
            "benchmarks/{sample}/metics.bwa.benchmark.txt"
    params:
        check=path.join(PATH_QC, "ReadStats", "{sample}.reads.q37"),       #these files are needed for the script
        check2=path.join(PATH_QC, "ReadStats", "{sample}.reads.aligned")   #these files are needed for the script 
    script:
        path.join(PATH_PIPELINE, 'scripts', 'WES_Combine_metrics.R')

rule Collect_Metrics:
    input:
        expand(path.join(PATH_QC, "Combined/{sample}_Combined_metrics.txt"), sample = AllFiles.keys())
    output:
        report(path.join(PATH_QC, "Combined_Metrics.txt"),category="Metrics")
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

#NEW ADDED
rule targetRegionFilter:
    input:
        bam=path.join(PATH_BAM, "{sample}_coordsorted.bam")
    output:
        exonbam=path.join(PATH_BAM, "{sample}_target_exons.bam")
    threads: config['all']['THREADS']
    conda:
        "envs/samtools.yaml"   #this env is based of  "envs/samtools2.yaml"
    params:
        targetbed=config["all"]["targets"]
    shell:
        "samtools view -@ {threads} -L {params.targetbed} {input.bam} -b -o {output.exonbam} &&"
        "samtools index {output.exonbam}"

#NEW ADDED
rule LoFreq_bam_indelQ:
    input:
        exonbam=path.join(PATH_BAM, "{sample}_target_exons.bam")
    output:
        indelQbam=temp(path.join(PATH_BAM, "{sample}_indelq_tmp.bam")),
        indelQbamSorted=path.join(PATH_BAM, "{sample}_indelq.bam")
    params:
       ref=config["all"]["REF"],
    threads: config['all']['THREADS']   
    conda:
        path.join(PATH_PIPELINE, "envs", "old_lofreq.yaml")    #this env is based of  "new_lofreq.yaml"
    log: path.join(PATH_LOG, ".LoFreq/{sample}_lofreq_indelq.txt")
    shell:
          "lofreq indelqual --dindel -f {params.ref} {input.exonbam} -o {output.indelQbam} 2> {log} &&"
          "samtools sort -o {output.indelQbamSorted} {output.indelQbam} 2>> {log}  &&"
          "samtools index {output.indelQbamSorted}"


rule LoFreq_paired:
    input:
        tumor=path.join(PATH_BAM, "{sample_paired}_indelq.bam"),
        normal = lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample_paired] + '_indelq.bam'),
    output:
        snps=path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_final.snvs.vcf.gz"),
        indels=path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_final.indels.vcf.gz"),
        check1=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_normal_relaxed.vcf.gz")),  		 #Is added extra, for rerun use
        check2=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_normal_relaxed.log")),           	 #Is added extra, for rerun use
        check3=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_normal_stringent.snvs.vcf.gz")),       #Is added extra, for rerun use
        check4=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_normal_stringent.indels.vcf.gz")),     #Is added extra, for rerun use
        check5=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_tumor_relaxed.vcf.gz")),         	 #Is added extra, for rerun use
        check6=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_tumor_relaxed.log")),        		 #Is added extra, for rerun use
        check7=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_tumor_stringent.snvs.vcf.gz")),        #Is added extra, for rerun use
        check8=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_tumor_stringent.indels.vcf.gz")),      #Is added extra, for rerun use
        check9=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_raw.snvs.vcf.gz")),         	 #Is added extra, for rerun use
        check10=temp(path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_raw.indels.vcf.gz")),         #Is added extra, for rerun use
    params:
        ref=config["all"]["REF"],
        threads=config["all"]["THREADS"],
        min_cov=config["LoFreq"]["min_cov"],
        min_mq=config["LoFreq"]["min_mq"],
        min_bq=config["LoFreq"]["min_bq"],
        min_alt_bq=config["LoFreq"]["min_alt_bq"],
        max_depth=config["LoFreq"]["max_depth"],
        sig=config["LoFreq"]["sig"],
        prefix=path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_")
    conda:
        path.join(PATH_PIPELINE, "envs", "lofreq.yaml")
    threads: config["all"]["THREADS"]
    log: path.join(PATH_LOG, "LoFreq", "{sample_paired}_lofreq.txt")
    benchmark:
            "benchmarks/{sample_paired}/Lofreq_paired.bwa.benchmark.txt"
    resources:    
          mem_mb=15000, 		#lambda wildcards, input, attempt: (input.size//1000000) * attempt * 6,
          disk_mb=15000, 		#lambda wildcards, input, attempt: (input.size//1000000) * attempt * 6
    shell:
        """
	lofreq somatic \
          -n {input.normal} -t {input.tumor} -f {params.ref} -o {params.prefix} \
        --min-cov {params.min_cov} \
        --threads {params.threads} --call-indels --verbose \
        2> {log}
        """

rule LoFreq_tumorOnly:
    input:
        tumor=path.join(PATH_BAM, "{sample_tumor}_indelq.bam")
    output:
        raw_snps=path.join(PATH_VAR, "LoFreq_tumorOnly", "{sample_tumor}_somatic_final.combined.vcf") 
    params:
        ref=config["all"]["REF"],
        min_cov=config["LoFreq"]["min_cov"],
        min_mq=config["LoFreq"]["min_mq"],
        min_bq=config["LoFreq"]["min_bq"],
        min_alt_bq=config["LoFreq"]["min_alt_bq"],
        max_depth=config["LoFreq"]["max_depth"],
        sig=config["LoFreq"]["sig"],
        prefix=path.join(PATH_VAR, "LoFreq_tumorOnly", "{sample_tumor}_")
    conda:
        path.join(PATH_PIPELINE, "envs", "lofreq.yaml")   #based on en "lofreq2.yaml"
    threads: config["all"]["THREADS"]
    benchmark:
            "benchmarks/{sample_tumor}/lofreq_tumorOnly.bwa.benchmark.txt"
    log: path.join(PATH_LOG, "LoFreq", "{sample_tumor}_lofreq.txt")
    shell:
        """
        lofreq call-parallel --pp-threads {threads} --call-indels --verbose {input.tumor} -f {params.ref} -o {output.raw_snps} \
        --min-cov {params.min_cov}  --min-mq {params.min_mq} --min-bq {params.min_bq}  --min-alt-bq {params.min_alt_bq}  --max-depth {params.max_depth}  --sig {params.sig}         2> {log}
        """

rule LoFreq_combine:
    input:
        snps=path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_final.snvs.vcf.gz"), 
        indels=path.join(PATH_VAR, "LoFreq_paired", "{sample_paired}_somatic_final.indels.vcf.gz")
    output:
        tmp=temp(path.join(PATH_VAR, "LoFreq", "{sample_paired}_somatic_final.combined.tmp.vcf")),
        out=path.join(PATH_VAR, "LoFreq", "{sample_paired}_somatic_final.combined.vcf")
    conda:
        "envs/lofreq.yaml"
    benchmark:
        "benchmarks/{sample_paired}/lofreq_combine.bwa.benchmark.txt"
    shell:
        """
	tabix -f {input.snps} &&    \
	tabix -f {input.indels} &&  \
        vcfcombine {input.snps} {input.indels} > {output.tmp} && \
        sed 's/\%//' {output.tmp} | sed 's/FREQ\,Number\=1\,Type\=String/FREQ\,Number\=1\,Type\=Float/' > {output.out}
	"""

def LoFreq_Mode(sample_L):
      if sample_L in Non_Paired:
              mode = "LoFreq_tumorOnly"
      elif sample_L in Tumor_Matched:
              mode = "LoFreq"
      return mode

rule LoFreq_readStatFilter:
    input:
        raw_vcf=lambda wildcards: expand(path.join(PATH_VAR, "{mode}", "{sample}_somatic_final.combined.vcf"),mode=LoFreq_Mode(wildcards.sample), sample=wildcards.sample )
    output:
        tmp_vcf=temp(path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_tmp.vcf")),
        unsorted_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_unsorted.vcf"),
        filtered_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_filt.vcf"),
    params:
        af_min=config["LoFreq"]["af_min"],
        cov_min=config["LoFreq"]["cov_min"],
        sb_alpha=config["LoFreq"]["sb_alpha"],
        hg19_dict=config["all"]["HG19_DICT"],
    conda:
        "envs/lofreq.yaml"
    log: path.join(PATH_LOG, "LoFreq/{sample}_lofreq_readStatFilter.txt")
    benchmark:
            "benchmarks/{sample}/lofreq_readStatFilter.bwa.benchmark.txt"
    shell:
        """
        lofreq filter --verbose --af-min {params.af_min} --cov-min {params.cov_min} --sb-alpha {params.sb_alpha} --sb-incl-indels -i {input.raw_vcf} -o {output.tmp_vcf}
        SnpSift filter -f {output.tmp_vcf} "(DP4[2]>2) & (DP4[3]>2) & ((na HRUN) | (HRUN<8))" > {output.unsorted_vcf}
        picard SortVcf I={output.unsorted_vcf} O={output.filtered_vcf} SD={params.hg19_dict} 2> {log}
        """

rule LoFreq_BlacklistFilter:
	input:
		filt_vcf=path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_filt.vcf") 
	output:
		blacklisted_vcf=path.join(PATH_VAR, "LoFreq/vcf/blacklisted/{sample}_blacklisted.vcf"),
		blacklisted_vcf_gz=path.join(PATH_VAR, "LoFreq/vcf/blacklisted/{sample}_blacklisted.vcf.gz"),
		blacklisted_csv=path.join(PATH_VAR, "LoFreq/vcf/blacklisted/{sample}_blacklisted.csv"),
		not_blacklisted_vcf=path.join(PATH_VAR, "LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf"),
		not_blacklisted_vcf_gz=path.join(PATH_VAR, "LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz"),
		not_blacklisted_csv=path.join(PATH_VAR, "LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.csv"),
	params:
		BED_blacklist=config["LoFreq"]["BED_blacklist"],
		Gene_blacklist=config["LoFreq"]["Gene_blacklist"],
		fields='CHROM POS REF ALT DP AF',
	conda:
		"envs/SnpSift.yaml"
	log:    path.join(PATH_LOG, "LoFreq/{sample}_lofreq_readStatFilter.txt") 	# "../logs/LoFreq/{sample}_lofreq_readStatFilter.txt" 
	shell:
		"""
		SnpSift intervals -i {input.filt_vcf} {params.BED_blacklist} {params.Gene_blacklist} > {output.blacklisted_vcf} 2> {log} &&
		SnpSift intervals -i {input.filt_vcf} -x {params.BED_blacklist} {params.Gene_blacklist} > {output.not_blacklisted_vcf} && 
		SnpSift extractFields -e "." {output.blacklisted_vcf} {params.fields} > {output.blacklisted_csv} &&
		SnpSift extractFields -e "." {output.not_blacklisted_vcf} {params.fields} > {output.not_blacklisted_csv} &&
		pbgzip -c {output.blacklisted_vcf} > {output.blacklisted_vcf_gz} &&
		tabix -s1 -b2 -e2 {output.blacklisted_vcf_gz} &&
		pbgzip -c {output.not_blacklisted_vcf} > {output.not_blacklisted_vcf_gz} &&
		tabix -s1 -b2 -e2 {output.not_blacklisted_vcf_gz} 
		"""

rule zip_vcf:
    input:
        path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf')
    output:
        vcf=path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf.gz'),
        idx=path.join(PATH_VAR, 'LoFreq/filtered/{sample}.vcf.gz.tbi')
    conda:
        "envs/lofreq.yaml"
    benchmark:
            "benchmarks/{sample}/zip_vcf.bwa.benchmark.txt"
    shell:
        """
        bgzip -c {input} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

#THIS RULE IS ONLY NEEDED IN SPECIAL CASE
#SPECIAL CASE IS: IF BAM HEADER IN FILES ARE WRONG
#RULES NEEDED ARE: "change_bam_header"  and  "Sambamba_sort2_bam_header"
rule change_bam_header:
    input: 
        tumor=path.join(PATH_BAM, "{sample_paired}_coordsorted.bam"),
        normal_bai=  lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample_paired] + '_coordsorted.bam.bai'),
        normal = lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample_paired] + '_coordsorted.bam'),
    output:
        bam=path.join(PATH_BAM, "{sample_paired}_normal_coordsorted.bam"),
        sam_header = temp(path.join(PATH_BAM, "{sample_paired}_new_header.sam")),
        sam_corrected = temp(path.join(PATH_BAM, "{sample_paired}_new_header_coordsorted.sam")),
    params:
        normal_name=getNormalSample,
        complete_name= lambda wildcards: pairs[wildcards.sample_paired]
    conda:
        "envs/samtools2.yaml"
    log: path.join(PATH_LOG, "sambamba/{sample_paired}_change_bam_header.log"),
    benchmark:
            "benchmarks/{sample_paired}/change_bam_header.bwa.benchmark.txt"
    shell:
        """
        samtools view -H {input.normal} > {output.sam_header}
        sed "s/SM:{params.normal_name}/SM:{params.complete_name}/" {output.sam_header} > {output.sam_corrected}
        samtools reheader {output.sam_corrected}  {input.normal} > {output.bam} 
        """

#THIS RULE IS ONLY NEEDED IN SPECIAL CASE
#SPECIAL CASE IS: IF BAM HEADER IN FILES ARE WRONG
#RULES NEEDED ARE: "change_bam_header"  and  "Sambamba_sort2_bam_header"
#IN RULE Mutect2_paired CHANGE INPUT STATEMENT
rule Sambamba_sort2_bam_header:
    input:
       dedup=path.join(PATH_BAM, "{sample_paired}_normal_coordsorted.bam"),
    output:
       coordsorted=path.join(PATH_BAM, "{sample_paired}_normal2_coordsorted.bam"),
       sortbai=path.join(PATH_BAM, "{sample_paired}_normal2_coordsorted.bam.bai"),
    params:
       tmpdir=PATH_TEMP,
    threads: config['all']['THREADS'],
    benchmark:
            "benchmarks/{sample_paired}/sambamba_sort2_bam_header.bwa.benchmark.txt"
    conda:
       "envs/picard.yaml"
    log: path.join(PATH_LOG, "sambamba/{sample_paired}_sort2.log"),
    shell:
       "sambamba sort -o {output.coordsorted} --tmpdir {params.tmpdir} -t {threads} -m 5000000000 {input.dedup} 2> {log} && "
       "sambamba index -t {threads} {output.coordsorted} {output.sortbai} 2>> {log}"
       

#THIS NEEDS TO BE RE DONE CAUSE TUMOR BAM IS NOW CHANGED BUT NORMAL IS NOT
#CHange input normal  depending on if you called rule Sambamba_sort2_bam_header
#INPUT:normal=path.join(PATH_BAM, "{sample_paired}_normal2_coordsorted.bam"),   needed if rules "change_bam_header"  and  "Sambamba_sort2_bam_header" are callded
#INPUT:normal=lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample_paired] + '_coordsorted.bam')  needed if bam files were not altered by  rules "change_bam_header"  and  "Sambamba_sort2_bam_header"
rule Mutect2_paired:
    input:
        tumor=path.join(PATH_BAM, "{sample_paired}_target_exons.bam"),  					# old bam file is path.join(PATH_BAM, "{sample_paired}_coordsorted.bam"),
        normal=lambda wildcards: path.join(PATH_BAM, pairs[wildcards.sample_paired] + '_coordsorted.bam'), 	# path.join(PATH_BAM, "{sample_paired}_normal2_coordsorted.bam"),     
    output:
        raw=path.join(PATH_VAR, "Mutect2_paired/{sample_paired}_Mutect2.vcf"),
        f1r2=path.join(PATH_VAR, "Mutect2_paired/{sample_paired}_f1r2.tar.gz")
    params:
        ref=config["all"]["REF"],
        pon=config["Mutect2"]["pon"],
        gnomad=config["Mutect2"]["gnomad"],
        interval=config["Mutect2"]["interval"],
        normal_name= lambda wildcards: pairs[wildcards.sample_paired], 
        tumor_name=getTumorSample                             							#This calls a function getTumorSample(wildcards)
    conda: "envs/gatk4.yaml"
    threads: config["all"]["THREADS"]
    benchmark:
            "benchmarks/{sample_paired}/mutect2_paired.bwa.benchmark.txt"
    resources:
            mem_mb=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8, 			#Can be altered if input size is smaller or larger	
            disk_mb=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8                       #Can be altered if input size is smaller or larger
    log:	path.join(PATH_LOG, "Mutect2/{sample_paired}_mutect2.txt")
    shell:
        """
        gatk Mutect2 --java-options "-Xmx16G"  -R {params.ref} -I {input.normal} -I {input.tumor} \
        -O {output.raw}  --tumor-sample {params.tumor_name}\
        --germline-resource {params.gnomad} \
        -pon {params.pon} -L {params.interval} --f1r2-tar-gz {output.f1r2} --max-mnp-distance 0 --native-pair-hmm-threads {threads} 2> {log}
        """


rule Mutect2_tumorOnly:
    input:
        tumor=path.join(PATH_BAM, "{sample_tumor}_target_exons.bam")  						#path.join(PATH_BAM, "{sample_tumor}_coordsorted.bam")
    output:
        raw=path.join(PATH_VAR, "Mutect2_tumorOnly/{sample_tumor}_Mutect2.vcf"),
        f1r2=path.join(PATH_VAR, "Mutect2_tumorOnly/{sample_tumor}_f1r2.tar.gz")
    params:
        ref=config["all"]["REF"],
        pon=config["Mutect2"]["pon"],
        gnomad=config["Mutect2"]["gnomad"],
        interval=config["Mutect2"]["interval"]
    conda: "envs/gatk4.yaml"
    threads: config["all"]["THREADS"]
    resources:
        mem_mb=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8,  				#Can be altered if input size is smaller or larger
        disk_mb=lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8				#Can be altered if input size is smaller or larger
    benchmark:
        "benchmarks/{sample_tumor}/mutect2_tumorOnly.bwa.benchmark.txt"
    log: path.join(PATH_LOG, "Mutect2/{sample_tumor}_mutect2.txt")
    shell:
        """
        gatk Mutect2 --java-options "-Xmx16G"   -R {params.ref} -I {input.tumor} \
        -O {output.raw} \
        --germline-resource {params.gnomad} \
        -pon {params.pon} -L {params.interval} --f1r2-tar-gz {output.f1r2} --max-mnp-distance 0 --native-pair-hmm-threads {threads} 2> {log}
        """


def get_mode(sample_f):
      if sample_f in Non_Paired:
           mode = "tumorOnly"
      elif sample_f in Tumor_Matched:
           mode = "paired"
      return mode

rule LearnReadOrientationModel:
    input:
        f1r2=lambda wildcards: expand(path.join(PATH_VAR, "Mutect2_{mode}/{sample}_f1r2.tar.gz"),mode=get_mode(wildcards.sample),sample=(wildcards.sample))
    output:
        obmodel=path.join(PATH_VAR, "Mutect2/{sample}_read-orientation-model.tar.gz")
    conda:
        "envs/gatk4.yaml"
    log:
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2_orientation.txt")
    benchmark:
            "benchmarks/{sample}/learnReadOrientationModel.bwa.benchmark.txt"
    resources:
        mem_mb=20000, 			#Can be altered if input size is smaller or larger   	 #lambda wildcards, input, attempt: (input.size//1000000) * attempt * 3,
        disk_mb=20000			#Can be altered if input size is smaller or larger	 #lambda wildcards, input, attempt: (input.size//1000000) * attempt * 3
    shell:
        """
        gatk --java-options "-Xmx4G"  LearnReadOrientationModel -I {input.f1r2} -O {output.obmodel} 2> {log}
	"""

# This rule is to determine contamination
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
    benchmark:
            "benchmarks/{sample}/calculateContamination.bwa.benchmark.txt"
    log: 
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2.txt")
    shell:
        """
        gatk GetPileupSummaries -I {input.exonbam} -O {output.pileup} -V {params.variants} -L {params.interval} 2>> {log} 
        gatk CalculateContamination -I {output.pileup} -tumor-segmentation {output.segments} -O {output.contamination} 2>> {log}
        """


#pass the learned read orientation model and contaminationtable in FilterMutectCalls
rule FilterMutectCalls:
    input:
        raw=lambda wildcards:expand(path.join(PATH_VAR, "Mutect2_{mode}/{sample}_Mutect2.vcf"), mode=get_mode(wildcards.sample), sample=wildcards.sample),
	obmodel=path.join(PATH_VAR, "Mutect2/{sample}_read-orientation-model.tar.gz"),
        segments=path.join(PATH_VAR, "Mutect2/{sample}_segments.table"),
        contamination=path.join(PATH_VAR, "Mutect2/{sample}_calculatecontamination.table"),
    output:
        filtered=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_filt.vcf"),
    params:
        ref=config["all"]["REF"],
        events=config["Mutect2"]["clustered"],      		#this is for the clustered events, default =2 but than artifacts cause that real events are filtered and SHM is filtered out
        af_min=config["Mutect2"]["af_min"],
        reads= config["Mutect2"]["reads_per_strand"],	
        min_reads= config["Mutect2"]["min_reads"],	
    conda:
        "envs/gatk4.yaml"
    log: 
        path.join(PATH_LOG, "Mutect2/{sample}_mutect2.txt")
    benchmark:
            "benchmarks/{sample}/FilterMutectCalls.bwa.benchmark.txt"
    shell:
        """
        gatk FilterMutectCalls -V {input.raw} -R {params.ref} \
        --max-events-in-region {params.events} --min-allele-fraction {params.af_min} \
        --tumor-segmentation {input.segments} --contamination-table {input.contamination} \
        --min-reads-per-strand {params.reads} --unique-alt-read-count {params.min_reads} \
        --ob-priors {input.obmodel} -O {output.filtered} 2>> {log} 
	"""

#changed snpSift filter to include germline or haplotyp if needed
#FILTER =~ 'PASS  is different from FILTER = 'PASS'
#SnpSift filter -f {input.filtered} "((FILTER =~ 'PASS' & (ROQ > 20 | TLOD > 20))" > {output.passed} &&
#SnpSift filter -f {input.filtered} "(FILTER = 'PASS' & (ROQ > 20 | TLOD > 20))" > {output.passed} &&
#SnpSift filter -f {input.filtered} "((FILTER = 'PASS' | FILTER =~ 'germline') & (ROQ > 20 | TLOD > 20))" > {output.passed} &&
#SnpSift filter -f {input.filtered} "((FILTER = 'PASS' | FILTER =~ 'germline' | FILTER =~ 'haplotype') & (ROQ > 20 | TLOD > 20))  > {output.passed} &&
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
    benchmark:
            "benchmarks/{sample}/mutect_passed.bwa.benchmark.txt"
    shell:
        """
        SnpSift filter -f {input.filtered} "((FILTER =~ 'PASS' | FILTER =~ 'germline') & (ROQ > 20 | TLOD > 20))" > {output.passed} &&
        SnpSift extractFields -e "." {output.passed} {params.fields}> {output.pass_csv} &&
        pbgzip -c {output.passed} > {output.passed_gz} &&
        tabix -s1 -b2 -e2 {output.passed_gz}
        """

rule Intersect_VariantCalls:
    input:
        lofreq_vcf=path.join(PATH_VAR, "LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf.gz"),
        Mutect2_vcf=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.vcf.gz"),
    output:
        intersect_vcf=path.join(PATH_VAR, "intersect/{sample}_intersect.vcf"),
        outersect_vcf=path.join(PATH_VAR, "intersect/{sample}_outersect.vcf")
    conda:
        "envs/SnpSift.yaml"
    log: path.join(PATH_LOG, "Intersect_variantCalls/{sample}_intersectCalls.txt")  
    benchmark:
            "benchmarks/{sample}/intersect_variantcalls.bwa.benchmark.txt"
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
    benchmark:
            "benchmarks/{sample}/annotate_variantcalls.bwa.benchmark.txt"
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
    conda:
        "envs/gatk4.yaml"
    log:	
        path.join(PATH_LOG, "Funcotator/{sample}_funcotator.txt")
    benchmark:
            "benchmarks/{sample}/funcotator.bwa.benchmark.txt"
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


def getTumor(wildcards):
    name=re.match('[a-zA-Z0-9\-]*', wildcards.sample).group(0)
    return(name)

rule vcf2maf:
    input:
        vcf = path.join(PATH_VAR, "funcotated/{sample}_funcotated.vcf"),
        fasta = config["vcf2maf"]["vep_fasta"],
        data_vep = config["vcf2maf"]["vep_data"],
    output:
        maf = path.join(PATH_VAR, "mafs/{sample}.maf"),
        check= temp(path.join(PATH_VAR, "funcotated/{sample}_funcotated.vep.vcf"))
    params:
        tumor=getTumor,
        file=path.join(PATH_VAR, "/funcotated/{sample}_funcotated.vep.vcf"),   
        new_path=path.join(PATH_VAR, "/funcotated/save_vep_vcf")
    conda: "envs/vcf2maf_new.yaml" 							#env based of vcf2_environment.yml"
    log:    
        path.join(PATH_LOG, "mafs/{sample}_vcf2maf.txt")
    threads: 10
    benchmark:
            "benchmarks/{sample}/vcf2maf.bwa.benchmark.txt"
    resources:
        mem_mb= lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8,
        disk_mb= lambda wildcards, input, attempt: (input.size//1000000) * attempt * 8
    shell:
        """
        mkdir -p variant/mafs
        z=$( echo $CONDA_PREFIX/bin )
        vcf2maf.pl --input-vcf {input.vcf} --output-maf {output.maf} --tumor-id {params.tumor} --ref-fasta {input.fasta} --vep-data {input.data_vep} --vep-path $z  --retain-info GT,AD,AF --vcf-tumor-id {params.tumor} 
        """ 

rule merge_mafs:
     input:
         #location_mafs = expand(path.join(PATH_VAR, 'mafs/')), # this is not used anymre
         check = expand(path.join(PATH_VAR, 'mafs/{sample}.maf'), sample=Tumor)   # This input is not used, But is stated here so that this rule runs after vcf2maf. 
     output:
         summaryPlot_name =  path.join(PATH_VAR, 'summaryPlot.pdf'),
         oncoplot_name = path.join(PATH_VAR, 'oncoplot.pdf'),
         check1 = path.join(PATH_VAR, 'merged.maf_clinicalData.txt'),
         check2 = path.join(PATH_VAR, 'merged.maf_geneSummary.txt'),
         check3 = path.join(PATH_VAR, 'merged.maf_maftools.maf'),
         check4 = path.join(PATH_VAR, 'merged.maf_sampleSummary.txt'),
         check5 = path.join(PATH_VAR, 'merged.maf_summary.txt'),
         #summary_name = path.join(PATH_VAR, 'merged.maf')
     log: 
         path.join(PATH_LOG, "merge_mafs/merge_mafs.txt") 
     conda: "envs/maftools.yaml"
     params:
         path_maf =  path.join(PATH_VAR, 'mafs'),
         summary_name = path.join(PATH_VAR, 'merged.maf'),
         #place = path.join(PATH_VAR, 'mafs', 'filter_maf_with_'),
     script:  
        path.join(PATH_PIPELINE,'scripts', 'maftools_oncoplot.R')


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
    benchmark:
            "benchmarks/{sample}/extractVcfFields_csv.bwa.benchmark.txt"
    shell:
        """
        SnpSift extractFields -e "." {input.funcotated} {params.fields} > {output.tmp_csv_snpsift} &&
        ./scripts/createFuncotationCsvFromFuncotatorVcf.sh {input.funcotated} > {output.tmp_csv_func} &&
        paste {output.tmp_csv_snpsift} {output.tmp_csv_func} > {output.funcotated_csv}
        """		


def cfm_directory(sample_cfm):
    if sample_cfm in Non_Paired:
           mode = "LoFreq_tumorOnly"
    elif sample_cfm in Tumor_Matched:
           mode = "LoFreq"
    return mode


rule mutation_counts:
   input:
      mutect_vcf_raw=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_filt.vcf"),
      mutect_vcf_filtered=path.join(PATH_VAR, "Mutect2/filtered/{sample}_Mutect2_passed.vcf"),
      lofreq_raw_vcf=lambda wildcards: expand(path.join(PATH_VAR, "{mode}/{sample_x}_somatic_final.combined.vcf"), mode=cfm_directory(wildcards.sample), sample_x=wildcards.sample),
      lofreq_filtered_vcf=path.join(PATH_VAR, "LoFreq/vcf/not_blacklisted/{sample}_not_blacklisted.vcf"), 			#path.join(PATH_VAR, "LoFreq/filtered/{sample}_lofreq_filt.vcf"),
      sample_intersect=path.join(PATH_VAR, "intersect/{sample}_intersect.vcf")
   output:
      file=path.join(PATH_VAR, "mutation_call", "{sample}_Found_mutations.txt")
   benchmark:
     "benchmarks/{sample}/mutation_counts.bwa.benchmark.txt"
   #params:
   #   new_dir =  path.join(PATH_VAR, "mutation_call")
   shell:
       """
       m_raw="$(cat {input.mutect_vcf_raw} | egrep -v  '^#' | wc -l)"
       m_fil="$(cat {input.mutect_vcf_filtered} | egrep -v  '^#' | wc -l)"
       l_raw="$(cat {input.lofreq_raw_vcf} | egrep -v  '^#' | wc -l )"
       l_fil="$(cat {input.lofreq_filtered_vcf} | egrep -v  '^#' | wc -l)"
       intersect="$(cat {input.sample_intersect} | egrep -v  '^#' | wc -l)"
       m_per="$(echo "scale=4;$m_fil/$m_raw*100" |bc)"
       l_per="$(echo "scale=4;$l_fil/$l_raw*100" |bc)" 
       m_inter="$(echo "scale=4;$intersect/$m_fil*100" |bc)"
       l_inter="$(echo "scale=4;$intersect/$l_fil*100" |bc)"
       
       mkdir -p variant/mutation_call 
 
       if (( $(echo "$m_inter > $l_inter" | bc -l) ))
       then
          Inter_percent=$m_inter
       else
          Inter_percent=$l_inter
       fi
       
       printf "%s\t%s\t%s\t%s\t%s\t%s\n%s_mutect:\t%s\t%s\t%.2f\t%s\t%.2f\n%s_lofreq:\t%s\t%s\t%.2f\t%s\t%.2f\n" "Sample_tool" "Raw" "Filtered" "(Raw/Filtered)%" "Intersect" "(Intersect/Filtered)%" "{wildcards.sample}" "$m_raw" "$m_fil" "$m_per" "$intersect" "$m_inter" "{wildcards.sample}" "$l_raw" "$l_fil" "$l_per" "$intersect" "$l_inter"  >  {output}
       printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Sample" "Raw_mutect2" "Filtered_Mutect2" "Raw_Lofreq" "Filtered_Lofreq" "Intersect" "Inter_percentage" >> {output}
       printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n" "{wildcards.sample}" "$m_raw" "$m_fil" "$l_raw" "$l_fil" "$intersect" "$Inter_percent" >> {output}

       """

rule merge_mutation_counts:
    input:
        expand(path.join(PATH_VAR, "mutation_call" ,"{sample}_Found_mutations.txt"), sample=Tumor)
    output:
         report(path.join(PATH_VAR, 'merged_mutation_counts.tsv'), category="Found mutations")
    shell:
        """
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "Sample" "Raw_mutect2" "Filtered_Mutect2" "Raw_Lofreq" "Filtered_Lofreq" "Intersect" "Inter_percentage" > {output}
        cat {input} |  egrep -v "_mutect" | egrep -v "_lofre" | egrep -v "Intersect"   >> {output}
        """

          
rule merge_benchmarks:
    output:
       output="benchmarks/merged/{sample}_merged_benchmark.tsv"
    params:
       path_file=["benchmarks/{sample}/"] ,
       full_path=["benchmarks/{sample}/*"] ,
       output= ["benchmarks/merged/{sample}_merged_benchmark.tsv"],
    shell:
        """
        bash scripts/script_merge_benchmark.sh "{params.path_file}" "{params.full_path}" "{output.output}" 
        """ 

