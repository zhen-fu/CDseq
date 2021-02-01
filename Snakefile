import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.15.0")

##### load config and sample sheets #####
configfile: "bin/config.yaml"

samples = pd.read_table("bin/samples.tsv").set_index("sample", drop=False)

##### target rules #####

rule all:
    input:
        # # mergeLanesAndRename_PE
        # expand("raw_data/{samples.sample}-R1.fastq.gz", samples=samples.itertuples()),
        # expand("raw_data/{samples.sample}-R2.fastq.gz", samples=samples.itertuples()),
        # # trim_galore_hardtrim5
        # expand("analysis/trim_reads/{samples.sample}-R1.50bp_5prime.fq.gz", samples=samples.itertuples()),
        # expand("analysis/trim_reads/{samples.sample}-R2.50bp_5prime.fq.gz", samples=samples.itertuples()),
        # # trim_galore
        # expand("analysis/trim_reads/{samples.sample}-R1.50bp_5prime_val_1.fq.gz", samples=samples.itertuples()),
        # expand("analysis/trim_reads/{samples.sample}-R2.50bp_5prime_val_2.fq.gz", samples=samples.itertuples()),
        # bwa_mem
        expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam", samples=samples.itertuples()),
        #expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam.bai", samples=samples.itertuples()),
        # extract_TLEN3
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam.bai", samples=samples.itertuples()),
       # mergeLanesAndRename_PE
        expand("raw_data/{samples.sample}-R1.fastq.gz", samples=samples.itertuples()),
        expand("raw_data/{samples.sample}-R2.fastq.gz", samples=samples.itertuples()),
        # trim_galore_hardtrim5
        expand("analysis/trim_reads/{samples.sample}-R1.80bp_5prime.fq.gz", samples=samples.itertuples()),
        expand("analysis/trim_reads/{samples.sample}-R2.80bp_5prime.fq.gz", samples=samples.itertuples()),
        # trim_galore
        expand("analysis/trim_reads/{samples.sample}-R1.80bp_5prime_val_1.fq.gz", samples=samples.itertuples()),
        expand("analysis/trim_reads/{samples.sample}-R2.80bp_5prime_val_2.fq.gz", samples=samples.itertuples()),
        # bwa_mem
        expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam", samples=samples.itertuples()),
        #expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam.bai", samples=samples.itertuples()),
        # extract_TLEN3
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam.bai", samples=samples.itertuples()),
        # extract_bed
       # "bin/2020-05-18_PFEG-20200427-CDseq.UVB-2.html"
        # deeptools_cov
        expand("analysis/deeptools/{samples.sample}.bw", samples=samples.itertuples()),
        # deeptools_heatmap
        expand("analysis/deeptools/{samples.sample}.compmat.gz", samples=samples.itertuples()),
        expand("analysis/deeptools/{samples.sample}.compmat.bed", samples=samples.itertuples()),
        expand("analysis/deeptools/{samples.sample}.profile.png", samples=samples.itertuples()),
        expand("analysis/deeptools/{samples.sample}.heatmap.png", samples=samples.itertuples()),

rule mergeLanesAndRename:
    input:
    output:
        "raw_data/{sample}-R1.fastq.gz",
        "raw_data/{sample}-R2.fastq.gz"
    log:
        "logs/mergeLanesAndRename/mergeLanesAndRename_PE-{sample}.log"
    benchmark:
        "benchmarks/mergeLanesAndRename/{sample}.bmk"
    resources:
        nodes = 1,
        threads = 1,
        mem_gb = 16
    envmodules:
        "bbc/R/R-3.6.0"
    script:
        "bin/mergeLanesAndRename.R"

rule trim_galore_hardtrim5:
    input:
        R1 = "raw_data/{sample}-R1.fastq.gz",
        R2 = "raw_data/{sample}-R2.fastq.gz",
    output:
        R1 = "analysis/trim_reads/{sample}-R1.80bp_5prime.fq.gz",
        R2 = "analysis/trim_reads/{sample}-R2.80bp_5prime.fq.gz",
    log:
        R1_o = "logs/trim_galore_hardtrim5/{sample}-R1.o",
        R1_e = "logs/trim_galore_hardtrim5/{sample}-R1.e",
        R2_o = "logs/trim_galore_hardtrim5/{sample}-R2.o",
        R2_e = "logs/trim_galore_hardtrim5/{sample}-R2.e",
    benchmark:
        "benchmarks/trim_galore_hardtrim5/{sample}.bmk",
    params:
        outdir="analysis/trim_reads/",
    envmodules:
        "bbc/trim_galore/trim_galore-0.6.0",
        "bbc/pigz/pigz-2.4",
    resources:
        nodes = 1,
        threads = 4,
        mem_gb=80,
    shell:
        """
        trim_galore \
        --hardtrim5 80 \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        {input.R1} 2> {log.R1_e} 1> {log.R1_o}

        trim_galore \
        --hardtrim5 80 \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        {input.R2} 2> {log.R2_e} 1> {log.R2_o}
        """
rule trim_galore:
    input:
        "analysis/trim_reads/{sample}-R1.80bp_5prime.fq.gz",
        "analysis/trim_reads/{sample}-R2.80bp_5prime.fq.gz",
    output:
        "analysis/trim_reads/{sample}-R1.80bp_5prime_val_1.fq.gz",
        "analysis/trim_reads/{sample}-R1.80bp_5prime_val_1_fastqc.html",
        "analysis/trim_reads/{sample}-R1.80bp_5prime_val_1_fastqc.zip",
        "analysis/trim_reads/{sample}-R1.80bp_5prime.fq.gz_trimming_report.txt",
        "analysis/trim_reads/{sample}-R2.80bp_5prime_val_2.fq.gz",
        "analysis/trim_reads/{sample}-R2.80bp_5prime_val_2_fastqc.html",
        "analysis/trim_reads/{sample}-R2.80bp_5prime_val_2_fastqc.zip",
        "analysis/trim_reads/{sample}-R2.80bp_5prime.fq.gz_trimming_report.txt",
    log:
        stdout="logs/trim_reads/{sample}.o",
        stderr="logs/trim_reads/{sample}.e"
    benchmark:
        "benchmarks/trim_reads/{sample}.bmk"
    params:
        outdir="analysis/trim_reads/",
    envmodules:
        "bbc/trim_galore/trim_galore-0.6.0",
        "bbc/pigz/pigz-2.4",
    resources:
        nodes = 1,
        threads = 4,
        mem_gb=80,
    shell:
        """
        trim_galore \
        --paired \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        -q 20 \
        --fastqc \
        {input} \
        2> {log.stderr} 1> {log.stdout}
        """

rule bwa_mem:
    input:
        R1 = "analysis/trim_reads/{sample}-R1.80bp_5prime_val_1.fq.gz",
        R2 = "analysis/trim_reads/{sample}-R2.80bp_5prime_val_2.fq.gz",
    output:
        sorted_dupremoved_q20_bam = "analysis/align/{sample}.sorted.dupremoved.q20.bam",
        #sorted_dupremoved_q20_bam_bai = "analysis/align/{sample}.sorted.dupremoved.q20.bam.bai"
    log:
        bwa =           "logs/align/bwa_mem.{sample}.log",
        samblaster =    "logs/align/samblaster.{sample}.log",
        samtools_view = "logs/align/samtools_view.{sample}.log",
        samtools_sort = "logs/align/samtools_sort.{sample}.log",
    benchmark:
        "benchmarks/bwa/{sample}.bmk"
    params:
        ref = config["ref"]["index"]
    envmodules: 
        "bbc/bwa/bwa-0.7.17",
        "bbc/samblaster/samblaster-0.1.24",
        "bbc/samtools/samtools-1.9",
    resources:
        nodes = 1,
        threads = 30,
        mem_gb = 300,
    shell:
        """
        bwa mem -t {resources.threads} {params.ref} {input.R1} {input.R2} 2> {log.bwa}| \
        samblaster --removeDups 2> {log.samblaster}| \
        samtools view -@ {resources.threads} -h -q 20 2> {log.samtools_view}| \
        samtools sort -@ {resources.threads} -O BAM -o {output.sorted_dupremoved_q20_bam} 2> {log.samtools_sort}
        """

rule TLEN3_bam:
    input:
        bam = "analysis/align/{sample}.sorted.dupremoved.q20.bam"
    output:
        header =         temp("analysis/extract_TLEN3/{sample}.TLEN3_header.sam"),
        TLEN3_noHeader = temp("analysis/extract_TLEN3/{sample}.TLEN3_noHeader.sam"),
        TLEN3_bam =           "analysis/extract_TLEN3/{sample}.TLEN3.bam",
        TLEN3_bai =           "analysis/extract_TLEN3/{sample}.TLEN3.bam.bai",
    log:
        extract_header = "logs/extract_TLEN3/{sample}.extract_header.log",
        extract_TLEN3 =  "logs/extract_TLEN3/{sample}.extract_TLEN3.log",
        return_header =  "logs/extract_TLEN3/{sample}.return_header.log",
        index =          "logs/extract_TLEN3/{sample}.index.log"
    benchmark:
        "benchmarks/extract_TLEN3/{sample}.bmk"
    resources:
        nodes = 1,
        threads = 16,
        mem_gb = 64
    envmodules:
        "bbc/samtools/samtools-1.9",
    shell:
        """
        samtools view -@ {resources.threads} -H {input.bam} 1> {output.header} 2> {log.extract_header}
        samtools view -@ {resources.threads} {input.bam} | awk '{{ if ($9 == 3 || $9 ==-3) {{ print }} }}' >> {output.TLEN3_noHeader} 2> {log.extract_TLEN3}
        cat {output.header} {output.TLEN3_noHeader} | samtools view -@ {resources.threads} -O BAM -o {output.TLEN3_bam} 2> {log.return_header}
        samtools index -@ {resources.threads} -b {output.TLEN3_bam} 2> {log.index}
        """

rule deeptools_cov:
    input:
        bam =             "analysis/extract_TLEN3/{sample}.TLEN3.bam",
    params:
        binsize =         50,
        norm_method =     "RPGC",
        eff_genome_size = 2407883318,
        ignore_chr =      "chrX chrY chrM",
    output:
        bigwig =          "analysis/deeptools/{sample}.bw",
    log:
        stdout =          "logs/deeptools_cov/{sample}.o",
        stderr =          "logs/deeptools_cov/{sample}.e",
    benchmark:
                          "benchmarks/deeptools/{sample}.txt",
    resources:
        nodes =           1,
        threads =         8,
        mem_gb =          100,
    envmodules:
                          "bbc/deeptools/deeptools-3.3.1",
    shell:
        # calulctate the coverage
        ## for SE data: we don't extend read length to fragment length because we'd have to make a up avalue for fragment length from the bioanalyzer data or something.
        ## for PE data: --extendReads used. "Reads with mates are always extended to match the fragment size defined by the two read mates.
        ##              Unmated reads, mate reads that map too far apart (>4x fragment length) or even map to different chromosomes are treated like single-end reads."
        ## eff. genome size was copied from https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html; we use the 'option 2(kmers)' method because bams were mapq filtered. Used the read length=50bp value
        """
        bamCoverage \
        -p {resources.threads} \
        --bam {input.bam} \
        -o {output.bigwig} \
        --binSize {params.binsize} \
        --normalizeUsing {params.norm_method} \
        --effectiveGenomeSize {params.eff_genome_size} \
        --ignoreDuplicates \
        --ignoreForNormalization {params.ignore_chr} \
        --extendReads \
        1>> {log.stdout} 2>> {log.stderr}
        """

rule deeptools_profile:
    input:
                     "analysis/deeptools/{sample}.bw",
    output:
        compmat =    "analysis/deeptools/{sample}.compmat.gz",
        compmatbed = "analysis/deeptools/{sample}.compmat.bed",
        profile =    "analysis/deeptools/{sample}.profile.png",
    log:             "logs/deeptools/{sample}.deeptools_profile.log",
    benchmark:
                     "benchmarks/deeptools/{sample}.deeptools_profile.txt",
    params:
        binSize =    50,
        after =      "3000",
        before =     "3000",
        body_len =   "5000",
        genes =      config["ref"]["annotation"],
        labels =     "{sample}",
    resources:
        nodes =      1,
        threads =    8,
        mem_gb =     100,
    envmodules:      "bbc/deeptools/deeptools-3.3.1",
    shell:
        """
        computeMatrix \
        scale-regions \
        -p {resources.threads} \
        --binSize {params.binSize} \
        -b {params.before} \
        -a {params.after} \
        --regionBodyLength {params.body_len} \
        -R {params.genes} \
        -S {input} \
        --samplesLabel {params.labels} \
        -o {output.compmat} \
        --skipZeros \
        --outFileSortedRegions {output.compmatbed} \
        2> {log}

        plotProfile \
        -m {output.compmat} \
        -out {output.profile} \
        2> {log}

        """

rule deeptools_heatmap:
    input:
        compmat =    "analysis/deeptools/{sample}.compmat.gz",
        compmatbed = "analysis/deeptools/{sample}.compmat.bed",
    output:
        heatmap =    "analysis/deeptools/{sample}.heatmap.png",
    log:             "logs/deeptools/{sample}.deeptools_heatmap.log",
    benchmark:
                     "benchmarks/deeptools/{sample}.deeptools_heatmap.txt",
    params:
        binSize =    50,
        after =      "3000",
        before =     "3000",
        body_len =   "5000",
        genes =      config["ref"]["annotation"],
        labels =     "{sample}",
    resources:
        nodes =      1,
        threads =    8,
        mem_gb =     100,
    envmodules:      "bbc/deeptools/deeptools-3.3.1",
    shell:
        """
        plotHeatmap \
        -m {input.compmat} \
        -out {output.heatmap} \
        2> {log}

        """


