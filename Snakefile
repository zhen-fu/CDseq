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
        # extract_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bed", samples=samples.itertuples()),
        # bed_sort_bgzip_tabix
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.sorted.bed.gz", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.sorted.bed.gz.tbi", samples=samples.itertuples()),
        # lenX_3prime_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.3prime.bed", len=[1,2], samples=samples.itertuples()),
        # lenX_5prime_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.5prime.bed", len=[1,2], samples=samples.itertuples()),
        # logo_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.logo.bed", samples=samples.itertuples()),
        # render_rmd
        expand("bin/{report}.html", report=config["report"]),

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
        #expand("analysis/align/{samples.sample}.sorted.dupremoved.q20.bam.bai", samples=samples.itertuples()),
        # extract_TLEN3
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bam.bai", samples=samples.itertuples()),
        # extract_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.bed", samples=samples.itertuples()),
        # bed_sort_bgzip_tabix
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.sorted.bed.gz", samples=samples.itertuples()),
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.sorted.bed.gz.tbi", samples=samples.itertuples()),
        # lenX_3prime_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.3prime.bed", len=[1,2], samples=samples.itertuples()),
        # lenX_5prime_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.5prime.bed", len=[1,2], samples=samples.itertuples()),
        # logo_bed
        expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.logo.bed", samples=samples.itertuples()),
        # render_rmd
        # "bin/2020-05-18_PFEG-20200427-CDseq.UVB-2.html"
        # deeptools_cov
        expand("analysis/deeptools/{samples.sample}.bw", samples=samples.itertuples()),
        # deeptools_heatmap
        expand("analysis/deeptools/{samples.sample}.compmat.gz", samples=samples.itertuples()),
        expand("analysis/deeptools/{samples.sample}.compmat.bed", samples=samples.itertuples()),
        expand("analysis/deeptools/{samples.sample}.profile.png", samples=samples.itertuples()),
        # "analysis/deeptools/compmat.gz",
        # "analysis/deeptools/compmat.bed",
        # "analysis/deeptools/heatmap.png",
        # getGCbigWig
        "analysis/getGCbigWig/mm10.GC.bw",
        # deeptools_profile_GC
        #"analysis/deeptools/profile.mm10.GC.png",
        # logo_bed
        "analysis/TES_motifs/mm10.TES.500bp.seq",
        "analysis/TES_motifs/mm10.TES.500bp.logo.bed",

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
        ref = "/primary/projects/bbc/references/mouse/indexes/mm10/bwa_gencode/mm10_gencode"
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

rule extract_TLEN3:
    input:
        TLEN3_bam = "analysis/extract_TLEN3/{sample}.TLEN3.bam",
    output:
        #sam =              temp("analysis/extract_TLEN3/{sample}.TLEN3.sam"),
        TLEN3_tmp =         temp("analysis/extract_TLEN3/{sample}.tmp"),
        seq =               temp("analysis/extract_TLEN3/{sample}.TLEN3.seq"),
        bed =                    "analysis/extract_TLEN3/{sample}.TLEN3.bed",
        sorted_bed_gz =          "analysis/extract_TLEN3/{sample}.TLEN3.sorted.bed.gz",
        sorted_bed_gz_tbi =      "analysis/extract_TLEN3/{sample}.TLEN3.sorted.bed.gz.tbi",

    params:
        ref =        expand("{ref}", ref=config["ref"]),
        tmp =        "analysis/extract_TLEN3/{sample}.TLEN3.tmp",
        sorted_bed = "analysis/extract_TLEN3/{sample}.TLEN3.sorted.bed",
    log:
        #get_sam =       "logs/extract_TLEN3/{sample}.get_sam.log",
        TLEN3_tmp =     "logs/extract_TLEN3/{sample}.TLEN3_tmp.log",
        get_bed =       "logs/extract_TLEN3/{sample}.get_bed.log",
        get_bed_name =  "logs/extract_TLEN3/{sample}.get_bed_name.log",
        get_seq =       "logs/extract_TLEN3/{sample}.get_seq.log",
        sort_bed =      "logs/extract_TLEN3/{sample}.sort_bed.log",

    benchmark:
        "benchmarks/extract_TLEN3/{sample}.bmk"
    resources:
        nodes = 1,
        threads = 12,
        mem_gb = 300,
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/bedtools/bedtools-2.29.2",
        "bbc/htslib/htslib-1.10.2",
    shell:
        """
        # make temp SAM for parsing; extract only TLEN=3 (TLEN=-3 are mates, i.e. duplicate entries)
        samtools view -@ {resources.threads} -h {input.TLEN3_bam} |\
        awk -v OFS='\t' '{{ if ($9 ==3) {{ print }} }}' 1> {output.TLEN3_tmp} 2> {log.TLEN3_tmp}

        # initiate BED [chr, start, stop]
        awk -v OFS='\t' '{{ print $3,$8-3,$8 }}' {output.TLEN3_tmp} 1> {output.bed} 2> {log.get_bed}

        # add name for row: [chr, start, stop, name]
        awk -v OFS='\t' '{{$(NF+1) = "cut_"i++}} 1' {output.bed} > {params.tmp} && mv {params.tmp} {output.bed}

        # get sequences from BED (5'-3' + strand by default)
        bedtools getfasta -tab -fi {params.ref} -bed {output.bed} -fo {output.seq} 2> {log.get_seq}

        # annotate strand based on middle nucleotide: T|C == (+)strand, A|G == (-)strand
        awk -v OFS='\t' '{{if( $2 ~ /.[TC]./ ) print $0,"0","+"; else print $0,"0","-"}}' {output.seq} 1> {params.tmp} && mv {params.tmp} {output.seq}

        # add seq, score, strand to BED: [chr, start, stop, name, (seq), score, strand]
        paste {output.bed} {output.seq} | cut -f 1,2,3,4,6,7,8 > {params.tmp} && mv {params.tmp} {output.bed}

        # merge name & seq fields: [chr, start, stop, name, score, strand]
        awk 'BEGIN {{FS=OFS="\t"}} {{ $4 = $4 "_" $5; print $1, $2, $3, $4, $6, $7 }}' {output.bed} > {params.tmp} && mv {params.tmp} {output.bed}

        # sort the BED file, bgzip it, index with tabix
        bedtools sort -i {output.bed} 1> {params.sorted_bed} 2> {log.sort_bed}
        bgzip -@ {threads} {params.sorted_bed}
        tabix {output.sorted_bed_gz}
        """

rule lengthX_3prime_bed:
    input:
        bed =            "analysis/extract_TLEN3/{sample}.TLEN3.bed",
    output:
        stop_plus = temp("analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.3prime.bed.stop_plus"),
        bed_plus =  temp("analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.3prime.bed.bed_plus"),
        bed =             "analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.3prime.bed",
    params:
        ref = expand("{ref}", ref=config["ref"]),
        tmp = "analysis/extract_TLEN3/{sample}.nt{num}.3prime.trinuc.tmp",
        nt =  "{num}"
    log:
        stop_plus =  "logs/lenX_3prime_bed/{sample}.nt{num}.3prime.edit_bed.log",
        get_seq =    "logs/lenX_3prime_bed/{sample}.nt{num}.3prime.get_seq.log",
        merge =      "logs/lenX_3prime_bed/{sample}.nt{num}.3prime.merge.log"
    resources:
        nodes = 1,
        threads = 1,
        mem_gb = 8,
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/bedtools/bedtools-2.29.2",
    shell:
        """
        # for (+) elements 'start'+1 and 'stop'+1, for (-) elements 'start'-1 and 'stop'-1.
        awk -v nt={params.nt} 'BEGIN {{FS=OFS="\t"}} \
        {{ if ($6 == "+") print $1, $2+1, $3+nt, $4, $5, $6; \
        else print $1, $2-nt, $3-1, $4, $5, $6 }}' {input.bed} 2> {log.stop_plus} 1> {output.stop_plus}

        # get sequences from BED in strand-specific manner.
        bedtools getfasta -s -tab -fi {params.ref} -bed {output.stop_plus} -fo {output.bed_plus} 2> {log.get_seq}

        # merge new coordinates [chr, start, stop] with strand-aware element sequence [name, score, strand]
        paste {output.stop_plus} {output.bed_plus} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{ $4 = "cut_"i++"_"$8; print $1, $2, $3, $4, $5, $6 }}' 1> {output.bed} 2> {log.merge}
        """

rule lengthX_5prime_bed:
    input:
        bed =            "analysis/extract_TLEN3/{sample}.TLEN3.bed",
    output:
        stop_plus = temp("analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.5prime.bed.stop_plus"),
        bed_plus =  temp("analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.5prime.bed.bed_plus"),
        bed =             "analysis/extract_TLEN3/{sample}.TLEN3.nt{num}.5prime.bed",
    params:
        ref = expand("{ref}", ref=config["ref"]),
        tmp = "analysis/extract_TLEN3/{sample}.nt{num}.5prime.trinuc.tmp",
        nt =  "{num}"
    log:
        stop_plus =  "logs/lenX_5prime_bed/{sample}.nt{num}.5prime.edit_bed.log",
        get_seq =    "logs/lenX_5prime_bed/{sample}.nt{num}.5prime.get_seq.log",
        merge =      "logs/lenX_5prime_bed/{sample}.nt{num}.5prime.merge.log"
    resources:
        nodes = 1,
        threads = 1,
        mem_gb = 8,
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/bedtools/bedtools-2.29.2",
    shell:
        """
        # for (+) elements 'start'-(nt-1) and 'stop', for (-) elements 'start'+(nt) and 'stop'.
        awk -v nt={params.nt} 'BEGIN {{FS=OFS="\t"}} \
        {{ if ($6 == "+") print $1, $2-(nt-1), $3, $4, $5, $6; \
        else print $1, $2, $3+(nt-1), $4, $5, $6 }}' {input.bed} 2> {log.stop_plus} 1> {output.stop_plus}

        # get sequences from BED in strand-specific manner.
        bedtools getfasta -s -tab -fi {params.ref} -bed {output.stop_plus} -fo {output.bed_plus} 2> {log.get_seq}

        # merge new coordinates [chr, start, stop] with strand-aware element sequence [name, score, strand]
        paste {output.stop_plus} {output.bed_plus} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{ $4 = "cut_"i++"_"$8; print $1, $2, $3, $4, $5, $6 }}' 1> {output.bed} 2> {log.merge}
        """

rule logo_bed:
    input:
        bed = "analysis/extract_TLEN3/{sample}.TLEN3.bed",
    output:
        bed_tmp =       temp("analysis/extract_TLEN3/{sample}.TLEN3.logo.bed_tmp"),
        logo_seq =                "analysis/extract_TLEN3/{sample}.TLEN3.logo.seq",
        logo_bed =           "analysis/extract_TLEN3/{sample}.TLEN3.logo.bed",
    params:
        ref = expand("{ref}", ref=config["ref"]),
    log:
        bed_tmp =  "logs/get_logo_seq/{sample}.bed_tmp.log",
        seq =      "logs/get_logo_seq/{sample}.seq_tmp.log",
        logo_bed = "logs/get_logo_seq/{sample}.logo.bed.log",
    resources:
        nodes = 1,
        threads = 1,
        mem_gb = 8,
    envmodules:
        "bbc/samtools/samtools-1.9",
        "bbc/bedtools/bedtools-2.29.2",
    shell:
        """
        # expand the start and stop coordinates to make a 20bp window around the cut site.
        awk 'BEGIN {{FS=OFS="\t"}} {{print $1, $2-9, $3+9, $4, $5, $6 }}' {input.bed} 1> {output.bed_tmp} 2> {log.bed_tmp}

        # get sequences from BED in strand-specific manner.
        bedtools getfasta -s -tab -fi {params.ref} -bed {output.bed_tmp} -fo {output.logo_seq} 2> {log.seq}

        # merge new coordinates [chr, start, stop] with strand-aware element sequence [name, score, strand]
        paste {output.bed_tmp} {output.logo_seq} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{ $4 = "cut_"i++"_"$8; print $1, $2, $3, $4, $5, $6 }}' 1> {output.logo_bed} 2> {log.logo_bed}
        """

rule render_rmd:
    input:
        nt_3prime = expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.3prime.bed", len=[1,2], samples=samples.itertuples()),
        nt_5prime = expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.nt{len}.5prime.bed", len=[1,2], samples=samples.itertuples()),
        logo_bed =  expand("analysis/extract_TLEN3/{samples.sample}.TLEN3.logo.bed", samples=samples.itertuples()),
        rmd =       expand("bin/{report}.Rmd", report=config["report"]),
    output:
        expand("bin/{report}.html", report=config["report"])
    params:
        ref = expand("{ref}", ref=config["ref"]),
    log:
        "logs/render_rmd.log",
    resources:
        nodes = "node095",
        threads = 1,
        mem_gb = 300,
    envmodules:
        # using the R on node095 and associated libraries required for rendering figures.
    shell:
        """
        # add rstudio-server pandoc installation to $PATH
        PATH=/usr/lib/rstudio-server/bin/pandoc:$PATH

        # render
        R -e 'rmarkdown::render("{input.rmd}")' 2> {log}
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

rule getGCbigWig:
    input:
    output:
        chrom_sizes =                   "analysis/getGCbigWig/mm10.chrom.sizes",
        mm10_gc5Base_wigVarStep_wig =   "analysis/getGCbigWig/mm10.gc5Base.wigVarStep.wig",
        GC_bw =                         "analysis/getGCbigWig/mm10.GC.bw",
    log:
        wget =                          "logs/getGCbigWig/wget.log",
        wigToBigWig =                   "logs/getGCbigWig/getGCbigWig.log",
    benchmark:                          "benchmarks/getGCbigWig/getGCbigWig.bmk",
    resources:
        nodes =                         1,
        threads =                       1,
        mem_gb =                        64,
    envmodules:                         "bbc/ucsc/ucsc-2020.06.11",
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes -O {output.chrom_sizes} 2> {log.wget}
        wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.gc5Base.wigVarStep.gz -O - | gunzip -c > {output.mm10_gc5Base_wigVarStep_wig} 2>> {log.wget}

        wigToBigWig {output.mm10_gc5Base_wigVarStep_wig} {output.chrom_sizes} {output.GC_bw} 2> {log.wigToBigWig}
        """

rule deeptools_computeMatrix_GC:
    input:
                     "analysis/getGCbigWig/mm10.GC.bw",
    output:
        compmat =    "analysis/deeptools/compmat.mm10.GC.gz",
        compmatbed = "analysis/deeptools/compmat.mm10.GC.bed",
    log:             "logs/deeptools/deeptools_computeMatrix.mm10.GC.log",
    benchmark:       "benchmarks/deeptools/deeptools_computeMatrix.mm10.GC.bmk",
    params:
        binSize =    50,
        after =      "3000",
        before =     "3000",
        body_len =   "5000",
        genes =      config["ref"]["annotation"],
        labels =     "mm10_%GC"
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
        """

rule TES_logo_bed:
    input:
        bed =        "analysis/TES_motifs/mm10.TES.500bp.bed",
    params:
        ref =        config["ref"]["sequence"],
    output:
        seq =        "analysis/TES_motifs/mm10.TES.500bp.seq",
        logo_bed =   "analysis/TES_motifs/mm10.TES.500bp.logo.bed",
    log:             "logs/TES_motifs/TES_motifs.log",
    benchmark:       "benchmarks/TES_motifs/TES_motifs.bmk",
    resources:
        nodes =      1,
        threads =    32,
        mem_gb =     100,
    envmodules:      "bbc/samtools/samtools-1.9",
                     "bbc/bedtools/bedtools-2.29.2",
    shell:
        """
        # get sequences from BED in strand-specific manner.
        bedtools getfasta -s -tab -fi {params.ref} -bed {input.bed} -fo {output.seq} 2> {log}

        # merge new coordinates [chr, start, stop] with strand-aware element sequence [name, score, strand]
        paste {input.bed} {output.seq} | \
        awk 'BEGIN {{FS=OFS="\t"}} {{ $4 = "cut_"i++"_"$8; print $1, $2, $3, $4, $5, $6 }}' 1> {output.logo_bed} 2>> {log}
        """
