import pandas as pd
from snakemake.utils import validate, min_version
##### set minimum snakemake version #####
min_version("5.15.0")

##### load config and sample sheets #####

samples = pd.read_table("src/samples.tsv").set_index("sample", drop=False)

##### target rules #####

rule all:
    input:
        # mergeLanesAndRename_PE
        expand("raw_data/{samples.sample}-R1.fastq.gz", samples=samples.itertuples()),
        expand("raw_data/{samples.sample}-R2.fastq.gz", samples=samples.itertuples()),
        # trim_galore_hardtrim5
        expand("analysis/trim_reads/{samples.sample}-R1.50bp_5prime.fq.gz", samples=samples.itertuples()),
        expand("analysis/trim_reads/{samples.sample}-R2.50bp_5prime.fq.gz", samples=samples.itertuples()),
        # trim_galore
        expand("analysis/trim_reads/{samples.sample}-R1.50bp_5prime_val_1.fq.gz", samples=samples.itertuples()),
        expand("analysis/trim_reads/{samples.sample}-R2.50bp_5prime_val_2.fq.gz", samples=samples.itertuples()),
        # bwa_mem
        expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam", samples=samples.itertuples()),
        # extract_TLEN3
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam.bai", samples=samples.itertuples()),
        # extract_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN_3.bed", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN_3_seq.txt", samples=samples.itertuples()),

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
        "src/mergeLanesAndRename.R"

rule trim_galore_hardtrim5:
    input:
        R1 = "raw_data/{sample}-R1.fastq.gz",
        R2 = "raw_data/{sample}-R2.fastq.gz",
    output:
        R1 = "analysis/trim_reads/{sample}-R1.50bp_5prime.fq.gz",
        R2 = "analysis/trim_reads/{sample}-R2.50bp_5prime.fq.gz",
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
        --hardtrim5 50 \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        {input.R1} 2> {log.R1_e} 1> {log.R1_o}

        trim_galore \
        --hardtrim5 50 \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        {input.R2} 2> {log.R2_e} 1> {log.R2_o}
        """
rule trim_galore:
    input:
        "analysis/trim_reads/{sample}-R1.50bp_5prime.fq.gz",
        "analysis/trim_reads/{sample}-R2.50bp_5prime.fq.gz",
    output:
        "analysis/trim_reads/{sample}-R1.50bp_5prime_val_1.fq.gz",
        "analysis/trim_reads/{sample}-R1.50bp_5prime_val_1_fastqc.html",
        "analysis/trim_reads/{sample}-R1.50bp_5prime_val_1_fastqc.zip",
        "analysis/trim_reads/{sample}-R1.50bp_5prime.fq.gz_trimming_report.txt",
        "analysis/trim_reads/{sample}-R2.50bp_5prime_val_2.fq.gz",
        "analysis/trim_reads/{sample}-R2.50bp_5prime_val_2_fastqc.html",
        "analysis/trim_reads/{sample}-R2.50bp_5prime_val_2_fastqc.zip",
        "analysis/trim_reads/{sample}-R2.50bp_5prime.fq.gz_trimming_report.txt",
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
        --three_prime_clip_R1 99 \
        --three_prime_clip_R2 99 \
        --output_dir {params.outdir} \
        --cores {resources.threads} \
        -q 20 \
        --fastqc \
        {input} \
        2> {log.stderr} 1> {log.stdout}
        """

rule bwa_mem:
    input:
        R1 = "analysis/trim_reads/{sample}-R1.50bp_5prime_val_1.fq.gz",
        R2 = "analysis/trim_reads/{sample}-R2.50bp_5prime_val_2.fq.gz",
    output:
        sorted_dupremoved_q20_bam = "analysis/align/{sample}.sorted.dupremoved.q20.bam"
    log:
        bwa =           "logs/align/bwa_mem.{sample}.log",
        samblaster =    "logs/align/samblaster.{sample}.log",
        samtools_view = "logs/align/samtools_view.{sample}.log",
        samtools_sort = "logs/align/samtools_sort.{sample}.log",
    benchmark:
        "benchmarks/bwa/{sample}.bmk"
    params:
        ref = "/primary/projects/bbc/references/mouse/indexes/mm10/bwa_gencode/mm10_gencode"
    envmodules:
        "bbc/bwa/bwa-0.7.17",
        "bbc/samblaster/samblaster-0.1.24",
        "bbc/samtools/samtools-1.9",
    resources:
        nodes = 1,
        threads = 8,
        mem_gb = 128,
    shell:
        """
        bwa mem -t {resources.threads} {params.ref} {input.R1} {input.R2} 2> {log.bwa}| \
        samblaster --removeDups 2> {log.samblaster}| \
        samtools view -@ {resources.threads} -h -q 20 2> {log.samtools_view}| \
        samtools sort -@ {resources.threads} -O BAM -o {output.sorted_dupremoved_q20_bam} 2> {log.samtools_sort}
        samtools index
        """

rule extract_TLEN3:
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

rule extract_TLEN3_bed:
    input:
        TLEN3_bam = "analysis/extract_TLEN3/{sample}.TLEN3.bam",
    output:
        sam =       temp("analysis/extract_TLEN3/{sample}.TLEN3.sam"),
        TLEN3_tmp = temp("analysis/extract_TLEN3/{sample}.tmp"),
        bed =            "analysis/extract_TLEN3/{sample}.TLEN_3.bed",
        seq =            "analysis/extract_TLEN3/{sample}.TLEN_3_seq.txt",
    params:
        ref = "/primary/projects/bbc/references/mouse/sequence/mm10/gdna/gencode/GRCm38.primary_assembly.genome.fa"
    log:
        get_sam =   "logs/extract_bed/{sample}.get_sam.log",
        TLEN3_tmp = "logs/extract_bed/{sample}.TLEN3_tmp.log",
        get_bed =   "logs/extract_bed/{sample}.get_bed.log",
        get_seq =   "logs/extract_bed/{sample}.get_seq.log",
    resources:
        nodes = 1,
        threads = 8,
        mem_gb = 32,
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/bedtools/bedtools-2.29.2",
    shell:
        """
        samtools view -@ {resources.threads} -h {input.TLEN3_bam} 1> {output.sam} 2> {log.get_sam}
        awk '{{ if ($9 ==3) {{ print }} }}' {output.sam} 1> {output.TLEN3_tmp} 2> {log.TLEN3_tmp}
        awk -v OFS='\t' '{{ print $3,$8-3,$8 }}' {output.TLEN3_tmp} 1> {output.bed} 2> {log.get_bed}
        bedtools getfasta -tab -s -fi {params.ref} -bed {output.bed} -fo {output.seq} 2> {log.get_seq}
        """
