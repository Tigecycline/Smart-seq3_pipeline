rule prepare_star_indices:
    input:
        ref_genome = ancient(config['reference']['genome']),
        gene_annotation = ancient(config['reference']['genes'])
    output:
        # use specified star index directory if available, otherwise default to a subdirectory next to genome fasta
        directory(join(dirname(config['reference']['genome']), 'star_index'))
    params:
        star_args = config['star_index_args']
    log: join(dirname(config['reference']['genome']), 'star_genome_generate.log')
    threads: workflow.cores
    conda: '../envs/star.yaml'
    shell:
        r'''
        STAR \
            {params.star_args} \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_genome} \
            --sjdbGTFfile {input.gene_annotation} \
            2> {log}
        '''


# rule link_fastq_files:
#     # so that umiextract produces the desired filename
#     input:
#         ancient(lambda wildcards: multiext(join(ilse_to_fastqdir[wildcards.seqid], '{fqid}/fastq/{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
#     output:
#         temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
#     shellz:
#         '''
#         ln -s {input[0]} {output[0]}
#         ln -s {input[1]} {output[1]}
#         '''


# rule merge_fastq:
#     input:
#         read1 = ancient(lambda wildcards: sample_to_read1[wildcards.sample]),
#         read2 = ancient(lambda wildcards: sample_to_read2[wildcards.sample])
#     output:
#         read1_merged = temp(join(config['outdir'], 'star_alignments/{sample}/{sample}_merged_R1.fastq.gz')),
#         read2_merged = temp(join(config['outdir'], 'star_alignments/{sample}/{sample}_merged_R2.fastq.gz'))
#     shell:
#         r'''
#         cat {input.read1} > {output.read1_merged}
#         cat {input.read2} > {output.read2_merged}
#         '''


rule extract_umi:
    input:
        read1 = lambda wildcards: join(fqid_to_dir[wildcards.fqid], '{fqid}/fastq/{fqid}_R1.fastq.gz'),
        read2 = lambda wildcards: join(fqid_to_dir[wildcards.fqid], '{fqid}/fastq/{fqid}_R2.fastq.gz')
    output:
        read1_umi = temp(join(config['outdir'], 'star_alignments/{sample}/{fqid}_R1_umiextract.fastq.gz')),
        read2_umi = temp(join(config['outdir'], 'star_alignments/{sample}/{fqid}_R2_umiextract.fastq.gz'))
        #temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1_umiextract.fastq.gz', '_R2_umiextract.fastq.gz'))
    params:
        outdir = join(config['outdir'], 'star_alignments/{sample}'),
        umiextract_args = config['umiextract_args']
    threads: min(4, workflow.cores)
    conda: '../envs/umite.yaml'
    shell:
        r'''
        umiextract \
            {params.umiextract_args} \
            -c {threads} \
            -1 {input.read1} \
            -2 {input.read2} \
            -d {params.outdir}
        '''


rule align_to_ref:
    input:
        read1 = lambda wildcards: expand(rules.extract_umi.output.read1_umi, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        read2 = lambda wildcards: expand(rules.extract_umi.output.read2_umi, sample=wildcards.sample, fqid=sample_to_fqid[wildcards.sample]),
        indices = ancient(rules.prepare_star_indices.output)
    output:
        join(config['outdir'], 'star_alignments/{sample}/{sample}_Aligned.out.bam')
    params:
        read1_comma = lambda wildcards, input: ','.join(input.read1),
        read2_comma = lambda wildcards, input: ','.join(input.read2),
        outprefix = join(config['outdir'], 'star_alignments/{sample}/{sample}_'),
        star_args = config['star_align_args']
    log: join(config['outdir'], 'star_alignments/{sample}/STAR_alignment.log')
    threads: min(4, workflow.cores)
    conda: '../envs/star.yaml'
    shell:
        r'''
        STAR \
            {params.star_args} \
            --runThreadN {threads} \
            --genomeDir {input.indices} \
            --genomeLoad LoadAndKeep \
            --readFilesIn {params.read1_comma} {params.read2_comma} \
            --outFileNamePrefix {params.outprefix} \
            > {log}
        '''


rule unload_star_genome:
    input:
        expand(rules.align_to_ref.output, sample=sample_to_fqid.keys())
    output:
        temp(touch(join(config['outdir'], 'star_genome_unloaded.flag')))
    conda: '../envs/star.yaml'
    params:
        indices = rules.prepare_star_indices.output
    shell:
        r'''
        STAR \
            --genomeLoad Remove \
            --genomeDir {params.indices}
        '''


rule sort_bam_by_query_name:
    input:
        rules.align_to_ref.output
    output:
        temp(join(config['outdir'], 'alignments/{sample}/{sample}_Aligned.qn_sorted.bam'))
    conda: '../envs/star.yaml'
    shell:
        'samtools cat {input} | samtools sort -n -o {output}' # umicount requires the BAM file to be sorted by query name (instead of genomic location)


rule parse_dump_GTF:
    input:
        ancient(config['reference']['genes'])
    output:
        join(dirname(config['reference']['genes']), 'umicount_GTF_dump.pkl')
    conda: '../envs/umite.yaml'
    shell:
        'umicount -g {input} --GTF_dump {output}'


rule count_umis:
    input:
        rules.unload_star_genome.output,
        bams = expand(rules.sort_bam_by_query_name.output, sample=sample_to_fqid.keys()),
        gtf_dump = ancient(rules.parse_dump_GTF.output)
    output:
        multiext(join(config['outdir'], 'umicount/umite'), '.D.tsv', '.RE.tsv', '.RI.tsv', '.UE.tsv', '.UI.tsv')
    params:
        outdir = join(config['outdir'], 'umicount'),
        umicount_args = config['umicount_args']
    log: join(config['outdir'], 'umicount/umicount.log')
    threads: min(4, workflow.cores)
    conda: '../envs/umite.yaml'
    shell:
        r'''
        umicount \
            {params.umicount_args} \
            -c {threads} \
            --GTF_skip_parse {input.gtf_dump} \
            --bams {input.bams} \
            -d {params.outdir} \
            -l {log}
        '''


# rule rename_umicount_files:
#     input:
#         rules.count_umis.output
#     output:
#         multiext(join(config['outdir'], f'umicount/{config['dataset']}_umite'), '.D.tsv', '.RE.tsv', '.RI.tsv', '.UE.tsv', '.UI.tsv')
#     params:
#         filename_prefix = config['dataset'],
#         samplename_suffix = '_Aligned.qn_sorted.bam'
#     shell:
#         'python3 scripts/rename_umicount_output.py --filename_prefix {params.filename_prefix} --samplename_suffix {params.samplename_suffix} {input}'


rule build_trsc_anndata:
    input:
        join(config['outdir'], 'umicount/umite.UE.tsv'),
        gtf_dump = ancient(rules.parse_dump_GTF.output)
    output:
        join(config['outdir'], f'{config['dataset']}.star_umite.h5ad')
    params:
        filename_prefix = config['dataset'],
        samplename_suffix = '_Aligned.qn_sorted.bam'
    conda: '../envs/anndata.yaml'
    shell:
        r'''
        python3 scripts/summarize_umicount_tsv.py \
            --filename_prefix {params.filename_prefix} \
            --samplename_suffix {params.samplename_suffix} \
            -g {input.gtf_dump} \
            {input[0]} \
            -o {output}
        '''