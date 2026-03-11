from snakemake.utils import min_version
min_version('9.0')


from os.path import basename, dirname, join
import yaml
import pandas as pd


configfile: 'snakeconfig.yaml'
with open('defaultconfig.yaml') as f:
    default_config = yaml.safe_load(f)
config = default_config | config


ilse_to_fastqdir = {}
sample_to_fastq = {}

for ilse_id, attr in config['ilse_info'].items():
    ilse_to_fastqdir[str(ilse_id)] = attr['fastqdir']

    if attr['metadata'].endswith('.xls') or attr['metadata'].endswith('.xlsx'):
        df_ilse = pd.read_excel(attr['metadata'], dtype=str)
    else: # if not xls, assume csv format (might need to be changed)
        df_ilse = pd.read_csv(attr['metadata'], dtype=str)
    for _, row in df_ilse.iterrows():
        if row['Sample Name'] in sample_to_fastq:
            sample_to_fastq[row['Sample Name']][str(ilse_id)] = row['Unique ID / Lane']
        else:
            sample_to_fastq[row['Sample Name']] = {str(ilse_id): row['Unique ID / Lane']}




wildcard_constraints:
    sn = r'\w+', # sample name
    seqid = r'[0-9]+', # sequencing ID
    fqid = r'[A-Z0-9-]+' # FASTQ ID








rule prepare_star_indices:
    input:
        ref_genome = ancient(config['ref_genome']['seq']),
        gene_annotation = ancient(config['ref_genome']['genes'])
    output:
        directory(join(dirname(config['ref_genome']['seq']), 'star_indices'))
    params:
        star_args = config['star_index_args']
    log: join(dirname(config['ref_genome']['seq']), 'star_genome_generate.log')
    threads: 24
    conda: 'envs/star.yaml'
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


# Running STAR with --genomeLoad LoadAndRemove does the same thing
# rule preload_genome:
#     input:
#         ancient(rules.prepare_star_indices.output)
#     output:
#         temp(touch('genome_loaded.flag'))
#     shell:
#         'STAR --genomeLoad LoadAndExit --genomeDir {input}'


rule link_fastq_files:
    # so that umiextract produces the desired filename
    input:
        ancient(lambda wildcards: multiext(join(ilse_to_fastqdir[wildcards.seqid], '{fqid}/fastq/{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
    output:
        temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1.fastq.gz', '_R2.fastq.gz'))
    shell:
        '''
        ln -s {input[0]} {output[0]}
        ln -s {input[1]} {output[1]}
        '''


rule extract_umi:
    input:
        rules.link_fastq_files.output
    output:
        temp(multiext(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}'), '_R1_umiextract.fastq.gz', '_R2_umiextract.fastq.gz'))
    params:
        outdir = join(config['outdir'], 'alignments/{sample}'),
        umiextract_args = config['umiextract_args']
    threads: 4
    conda: 'envs/umite.yaml'
    shell:
        r'''
        umiextract \
            {params.umiextract_args} \
            -c {threads} \
            -1 {input[0]} \
            -2 {input[1]} \
            -d {params.outdir}
        '''


rule align_to_ref:
    input:
        fq_with_umi = rules.extract_umi.output,
        indices = ancient(rules.prepare_star_indices.output)
    output:
        temp(join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}_Aligned.out.bam'))
    params:
        outprefix = join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}_'),
        star_args = config['star_align_args']
    log: join(config['outdir'], 'alignments/{sample}/{seqid}_{fqid}.STAR.log')
    threads: 4
    conda: 'envs/star.yaml'
    shell:
        r'''
        STAR \
            {params.star_args} \
            --runThreadN {threads} \
            --genomeDir {input.indices} \
            --readFilesIn {input.fq_with_umi} \
            --outFileNamePrefix {params.outprefix} \
            > {log}
        '''


rule merge_sort_aligned_reads:
    input:
        lambda wildcards: [join(config['outdir'], f'alignments/{{sample}}/{seqid}_{fqid}_Aligned.out.bam') for seqid, fqid in sample_to_fastq[wildcards.sample].items()]
    output:
        temp(join(config['outdir'], 'alignments/{sample}/{sample}_RNA.qn_sorted.bam'))
    conda: 'envs/star.yaml'
    shell:
        'samtools cat {input} | samtools sort -n -o {output}' # umicount requires the BAM file to be sorted by query name (instead of genomic location)


rule parse_dump_GTF:
    input:
        ancient(config['ref_genome']['genes'])
    output:
        join(dirname(config['ref_genome']['genes']), 'umicount_GTF_dump.pkl')
    conda: 'envs/umite.yaml'
    shell:
        'umicount -g {input} --GTF_dump {output}'


rule count_umi_rename_columns:
    input:
        bams = expand(join(config['outdir'], 'alignments/{sn}/{sn}_RNA.qn_sorted.bam'), sn=sample_to_fastq.keys()),
        gtf_dump = ancient(rules.parse_dump_GTF.output)
    output:
        multiext(join(config['outdir'], 'umicount/umite'), '.D.tsv', '.RE.tsv', '.RI.tsv', '.UE.tsv', '.UI.tsv')
    params:
        outdir = join(config['outdir'], 'umicount'),
        umicount_args = config['umicount_args']
    log: join(config['outdir'], 'umicount/umicount.log')
    threads: 4
    conda: 'envs/umite.yaml'
    shell:
        r'''
        umicount \
            {params.umicount_args} \
            -c {threads} \
            --GTF_skip_parse {input.gtf_dump} \
            --bams {input.bams} \
            -d {params.outdir} \
            2> {log}
        '''


rule rename_umicount_files:
    input:
        rules.count_umi_rename_columns.output
    output:
        multiext(join(config['outdir'], f'umicount/{config['dataset']}_umite'), '.D.tsv', '.RE.tsv', '.RI.tsv', '.UE.tsv', '.UI.tsv')
    params:
        filename_prefix = config['dataset'],
        samplename_suffix = '_RNA.qn_sorted.bam'
    shell:
        'python3 scripts/rename_umicount_output.py --filename_prefix {params.filename_prefix} --samplename_suffix {params.samplename_suffix} {input}'


rule build_trsc_anndata:
    input:
        join(config['outdir'], f'umicount/{config['dataset']}_umite.UE.tsv'),
        gene_annotation = ancient(config['ref_genome']['genes'])
    output:
        join(config['outdir'], f'{config['dataset']}_transcriptome.h5ad')
    conda: 'envs/anndata.yaml'
    shell:
        'python3 scripts/build_trsc_anndata.py -g {input.gene_annotation} {input[0]} -o {output}'


rule all:
    default_target: True
    input:
        rules.build_trsc_anndata.output
    output:
        join(config['outdir'], 'last_success_config.yaml')
    run:
        from datetime import datetime
        with open(output[0], 'w') as f:
            f.write(f'# These are configurations used for the most recent successful run of the pipeline, which was {datetime.today().isoformat(sep=' ', timespec='seconds')}.\n')
            yaml.dump(config, f, sort_keys=False)