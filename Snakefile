# Snakemake simple RNA-Seq pipeline

import os, sys, glob
import pandas as pd

sample_sheet = pd.read_table('/oak/stanford/groups/pritch/users/jweinstk/perturbation_data/metadata/sample_meta_data_2022_09_01.tsv', index_col = False)
rna_samples = sample_sheet.Sample.tolist()

fastq_dir = "/oak/stanford/groups/pritch/users/jweinstk/perturbation_data/fastq"
output_base_dir = "/oak/stanford/groups/pritch/users/jweinstk/perturbation_data/rnaseq_pipeline/output"

def create_directory_if_not_exist(dir)
    if not os.path.exists(dir):
        os.makedirs(dir)

    
dirs = [
    output_base_dir,
    os.path.join(output_base_dir, "qc"),
    os.path.join(output_base_dir, "qc", "rseqc"),
    os.path.join(output_base_dir, "qc", "fastqc"),
    os.path.join(output_base_dir, "qc", "preseq"),
    os.path.join(output_base_dir, "qc", "read_stats"),
    os.path.join(output_base_dir, "bam", "dedup"),
    os.path.join(output_base_dir, "bam", "star"),
    os.path.join(output_base_dir, "counts"),
    os.path.join(output_base_dir, "error_files"),
    os.path.join(output_base_dir, "biotype_final")
]

# all rule
rule all:
    input:
        "output/multiqc_report.html",
        "output/counts/dedup_counts.txt",
        "output/counts/all_counts.txt"

rule multiqc:
    input:
        expand("output/bam/dedup/{sample_label}.dedup.bam", sample_label = rna_samples),
        "output/counts/dedup_counts.txt",
        "output/counts/all_counts.txt",
        expand("output/qc/rseqc/{sample_label}_distribution.txt", sample_label = rna_samples),
        expand("output/qc/rseqc/{sample_label}_strandedness.txt", sample_label = rna_samples),
        expand("output/qc/rseqc/{sample_label}.GC.xls", sample_label = rna_samples),
        expand("output/qc/rseqc/{sample_label}.pos.DupRate.xls", sample_label = rna_samples),
        expand("output/qc/rseqc/{sample_label}.rawCount.xls", sample_label = rna_samples),
        expand("output/biotype_final/{sample_label}_unique_biotype_counts_mqc.txt", sample_label = rna_samples),
        expand("output/qc/fastqc/{sample_label}_umi_trimmed_fastqc.html", sample_label = rna_samples),
        expand("output/qc/preseq/{sample_label}.extrapolated_yield.txt", sample_label = rna_samples),
        expand("output/qc/read_stats/{sample_label}_deduplicated_mqc.csv", sample_label = rna_samples)

    output:
        "output/multiqc_report.html"
    params:
        error_out_file = "error_files/multiqc",
        run_time="1:00:00",
        cores="1",
        memory="4000",
        job_name="multiqc",
    benchmark:
        "benchmarks/rseqc/rseqc.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        "multiqc -f output/ error_files/ logs/ --ignore output/biotype/ -o output/"

# UMI-tools to extract UMI and add it to read header
# Also remove linker sequence
rule umiheader:
    input:
        # fastq = lambda wildcards: glob.glob('original_fastq/*/' + sample_sheet.loc[sample_sheet.Sample == wildcards.sample_label]["Davis_ID"].values[0] + '_UMI*.fastq.gz')
        fastq = lambda wildcards: sample_sheet.loc[sample_sheet.Sample == wildcards.sample_label]["path"].values[0]
    output:
        fastq = "output/fastq/umi/{sample_label}_umi.fastq.gz",
        log = "logs/umi_tools/{sample_label}_header.txt"
    params:
        error_out_file = "error_files/{sample_label}_umi_header",
        run_time="02:00:00",
        cores="1",
        memory="40000",
        job_name="umi_header"
    conda:
        "envs/rna_seq.yaml"
    benchmark: 
        "benchmarks/umi_tools/{sample_label}_header.txt"
    shell:
        "zcat {input.fastq} | umi_tools extract -S {output.fastq} --log={output.log} "
        "--extract-method=regex --bc-pattern='(?P<umi_1>.{{6}})(?P<discard_1>.{{4}}).*'"

# cutadapt
rule cutadapt:
    input:
        fastq = "output/fastq/umi/{sample_label}_umi.fastq.gz"
    output:
        fastq = "output/fastq/trimmed/{sample_label}_umi_trimmed.fastq.gz"
    params:
        error_out_file = "error_files/{sample_label}_cutadapt",
        run_time="2:00:00",
        cores="1",
        memory="6000",
        job_name="cutadapt"
    benchmark:
        "benchmarks/cutadapt/{sample_label}.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        """
        cutadapt -a "A{{30}}" -a "T{{30}}" -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --minimum-length 20 -o {output.fastq} {input.fastq}
        """

# map
rule star:
    input:
        "output/fastq/trimmed/{sample_label}_umi_trimmed.fastq.gz"
    output:
        bam = "output/bam/star/{sample_label}.unique.Aligned.sortedByCoord.out.bam",
        idx = "output/bam/star/{sample_label}.unique.Aligned.sortedByCoord.out.bam.bai"
    params:
        error_out_file = "error_files/{sample_label}_star",
        run_time="24:00:00",
        cores="1",
        memory="40000",
        job_name="star",
        prefix = "output/bam/star/{sample_label}.unique."
    benchmark:
        "benchmarks/star/{sample_label}.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        "STAR --genomeDir /oak/stanford/groups/pritch/users/jake/genome/human/star_human "
        "--readFilesIn {input} "
        "--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate "
        "--outFilterMultimapNmax 1 --outFileNamePrefix {params.prefix}; "
        "samtools index {output}"

# deduplicate based on UMIs
rule deduplicate:
    input:
        rules.star.output.bam
    output:
        bam = "output/bam/dedup/{sample_label}.dedup.bam",
        idx = "output/bam/dedup/{sample_label}.dedup.bam.bai",
        #log = "logs/umi_tools/{sample_label}_dedup_stats.txt_per_umi.tsv"
    params:
        error_out_file = "error_files/{sample_label}_umi_dedup",
        run_time="12:00:00",
        cores="1",
        memory="40000",
        job_name="umi_dedup",
        prefix = "logs/umi_tools/{sample_label}_dedup_stats"
    conda:
        "envs/rna_seq.yaml"
    benchmark: 
        "benchmarks/umi_tools/{sample_label}_dedup.txt"
    shell:
        "umi_tools dedup -I {input} --output-stats={params.prefix} -S {output.bam}; "
        "samtools index {output.bam}"

# featurecounts
rule featureCounts:
    input:
        dedup_bams = expand("output/bam/dedup/{sample_label}.dedup.bam", sample_label = rna_samples),
        all_bams = expand("output/bam/star/{sample_label}.unique.Aligned.sortedByCoord.out.bam", sample_label = rna_samples)
    output:
        dedup_counts = "output/counts/dedup_counts.txt",
        all_counts = "output/counts/all_counts.txt"
    params:
        error_out_file = "error_files/featurecounts",
        run_time="4:00:00",
        cores="1",
        memory="6000",
        job_name="featureCounts"
    benchmark:
        "benchmarks/featurecounts/featureCounts.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        "featureCounts -s 1 -a /oak/stanford/groups/pritch/users/jake/genome/human/gencode.v35.basic.annotation.gtf "
        "-o output/counts/dedup_counts.txt {input.dedup_bams}; "
        "featureCounts -s 1 -a /oak/stanford/groups/pritch/users/jake/genome/human/gencode.v35.basic.annotation.gtf "
        "-o output/counts/all_counts.txt {input.all_bams}; "

# biotype
rule biotype:
    input:
        rules.deduplicate.output.bam
    output:
        biotype_intermediate = "output/biotype/{sample_label}_biotype_unique.txt",
        biotype_final = "output/biotype_final/{sample_label}_unique_biotype_counts_mqc.txt"
    params:
        error_out_file = "error_files/{sample_label}_biotype",
        run_time="4:00:00",
        cores="1",
        memory="6000",
        job_name="biotype"
    benchmark:
        "benchmarks/biotype/{sample_label}.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        "featureCounts -s 1 -a /oak/stanford/groups/pritch/users/jake/genome/human/gencode.v35.basic.annotation.gtf "
        "-g gene_type -o {output.biotype_intermediate} {input}; "
        """echo "# section_name: 'featureCounts Biotype'" > {output.biotype_final}; """
        "cut -f 1,7 {output.biotype_intermediate} | tail -n +3 >> {output.biotype_final}"

# QC ---------------------------------

# fastqc
rule fastqc:
    input:
        rules.cutadapt.output.fastq
    output:
        "output/qc/fastqc/{sample_label}_umi_trimmed_fastqc.html",
        # stuff we don't really care about but want to eliminate when run is botched
        "output/qc/fastqc/{sample_label}_umi_trimmed_fastqc.zip"
    params:
        error_out_file = "error_files/{sample_label}_fastqc",
        run_time="00:15:00",
        cores="1",
        memory="6000",
        job_name="fastqc",
    conda:
        "envs/rna_seq.yaml"
    benchmark: 
        "benchmarks/fastqc/{sample_label}.txt"
    shell:
        "fastqc {input} --outdir=output/qc/fastqc/"

# Use Preseq to estimate library complexity
rule estimate_library_complexity:
    input:
        bam = rules.star.output.bam
    output:
        lc = "output/qc/preseq/{sample_label}.extrapolated_yield.txt"
    params:
        error_out_file = "error_files/{sample_label}_estimate_lc",
        run_time = "1:00:00",
        cores = "1",
        memory = "8000",
        job_name = "lc_extrap"
    conda:
        "envs/rna_seq.yaml"
    benchmark: "benchmarks/preseq/{sample_label}.txt"
    shell:
        "preseq lc_extrap -P -o {output.lc} -B {input.bam}" 

# Calculate how many reads were deduplicated
rule dedup_stats:
    input:
        unfiltered_bam = rules.star.output.bam,
        unfiltered_idx = rules.star.output.idx,
        dedup_bam = rules.deduplicate.output.bam,
        dedup_idx = rules.deduplicate.output.idx
    output:
        deduplicated = "output/qc/read_stats/{sample_label}_deduplicated_mqc.csv"
    params:
        error_out_file="error_files/{sample_label}_read_stats",
        run_time="0:20:00",
        cores="1",
        memory="6000",
        job_name="read_stats"
    threads: 1
    conda:
        "envs/rna_seq.yaml"
    shell:
        """
        echo "# parent_id: custom_section" > {output.deduplicated};
        echo "# parent_name: 'Deduplication stats'" >> {output.deduplicated};
        echo "# parent_description: 'Percent of reads remaning after deduplication based on UMI'" >> {output.deduplicated};
        echo "# id: custom_bargraph_csv" >> {output.deduplicated};
        echo "# section_name: 'Deduplication stats'" >> {output.deduplicated};
        echo "# description: 'Percent of reads remaning after deduplication based on UMI'" >> {output.deduplicated};
        echo "# format: 'csv'" >> {output.deduplicated};
        echo "# plot_type: 'bargraph'" >> {output.deduplicated};
        echo "# pconfig:" >> {output.deduplicated};
        echo "#    id: 'custom_bargraph_w_header'" >> {output.deduplicated};
        echo "#    title: Deduplication stats" >> {output.deduplicated};
        echo "#    ylab: '% reads remaining'" >> {output.deduplicated};

        mappedReads=$(samtools idxstats {input.unfiltered_bam} | awk '{{SUM += $3}} END {{print SUM}}');
        dedupReads=$(samtools idxstats {input.dedup_bam} | awk '{{SUM += $3}} END {{print SUM}}');
        dedupPercent=$(bc <<< "scale=2;$dedupReads/$mappedReads*100");
        echo "Post-deduplication,"${{dedupPercent}} >> {output.deduplicated};
        """

# rseqc
rule rseqc:
    input:
        bam = rules.deduplicate.output.bam,
        idx = rules.deduplicate.output.idx
    output:
        distribution = "output/qc/rseqc/{sample_label}_distribution.txt",
        strandedness = "output/qc/rseqc/{sample_label}_strandedness.txt",
        gc = "output/qc/rseqc/{sample_label}.GC.xls",
        duplication = "output/qc/rseqc/{sample_label}.pos.DupRate.xls",
        saturation = "output/qc/rseqc/{sample_label}.rawCount.xls"
    params:
        error_out_file = "error_files/{sample_label}_rseqc",
        run_time="24:00:00",
        cores="1",
        memory="12000",
        job_name="rseqc",
        prefix = "output/qc/rseqc/{sample_label}"
    benchmark:
        "benchmarks/rseqc/{sample_label}.txt"
    conda:
        "envs/rna_seq.yaml"
    shell:
        "geneBody_coverage.py -r /oak/stanford/groups/pritch/users/jake/genome/rseqc/hg38_Gencode_V28.bed "
        " -i {input.bam} -o {params.prefix}; "
        "read_distribution.py -i {input} -r /oak/stanford/groups/pritch/users/jake/genome/rseqc/hg38_Gencode_V28.bed > {output.distribution}; "
        "read_GC.py -i {input.bam} -o {params.prefix}; "
        "read_duplication.py -i {input.bam} -o {params.prefix}; "
        "infer_experiment.py -r /oak/stanford/groups/pritch/users/jake/genome/rseqc/hg38_Gencode_V28.bed -i {input.bam} > {output.strandedness}; "
        "RPKM_saturation.py -i {input.bam} -o {params.prefix} -r /oak/stanford/groups/pritch/users/jake/genome/rseqc/hg38_Gencode_V28.bed"

