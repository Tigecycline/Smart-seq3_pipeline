rule prepare_salmon_index:
    input:
        genome = ancient(config['reference']['genome']),
        transcriptome = ancient(config['reference']['transcriptome']),
    output:
        #decoys = temp(join(dirname(config['reference']['genome']), 'decoys.txt')),
        #gentrome = temp(join(dirname(config['reference']['genome']), 'gentrome.fa.gz')),
        index = directory(join(dirname(config['reference']['genome']), 'salmon_index')),
    params:
        salmon_index_args = config['salmon_index_args']
    log:
        join(dirname(config['reference']['genome']), 'salmon_index/salmon_index.log')
    threads: workflow.cores
    conda: '../envs/salmon.yaml'
    shell:
        # step 1: extract chromosome/contig from the reference genome and use as decoy
        # step 2: append decoy entries (which is the entire genome) to the transcriptome
        # step 3: feed the joint transcriptome+genome to salmon index
        # r'''
        # zcat -f {input.genome} | grep "^>" | cut -d " " -f 1 | tr -d ">" > {output.decoys}
        # cat {input.transcriptome} {input.genome} > {output.gentrome}
        # salmon index -p {threads} {params.salmon_index_args} -t {output.gentrome} -d {output.decoys} -i {output.index} &> {log}
        # '''
        r'salmon index -p {threads} {params.salmon_index_args} -t {input.transcriptome} -i {output.index} &> {log}'


rule salmon_quantify:
    input:
        read1 = ancient(lambda wildcards: [join(fqid_to_dir[fqid], f'{fqid}/fastq/{fqid}_R1.fastq.gz') for fqid in sample_to_fqid[wildcards.sample]]),
        read2 = ancient(lambda wildcards: [join(fqid_to_dir[fqid], f'{fqid}/fastq/{fqid}_R2.fastq.gz') for fqid in sample_to_fqid[wildcards.sample]]),
        index = ancient(rules.prepare_salmon_index.output.index)
    output:
        directory(join(config['outdir'], 'salmon_quants/{sample}'))
    params:
        salmon_quant_args = config['salmon_quant_args']
    log: join(config['outdir'], 'salmon_quants/{sample}/salmon_quantify.log')
    threads: min(4, workflow.cores)
    conda: '../envs/salmon.yaml'
    shell:
        r'salmon quant {params.salmon_quant_args} -p {threads} -i {input.index} -1 {input.read1} -2 {input.read2} -o {output} &> {log}'


rule summarize_quants:
    input:
        quant_dirs = expand(rules.salmon_quantify.output, sample=sample_to_fqid.keys()),
        gtf = ancient(config['reference']['genes'])
    output:
        join(config['outdir'], f'{config["dataset"]}.salmon.h5ad')
    conda: '../envs/anndata.yaml'
    shell:
        r'python3 scripts/summarize_salmon_quants.py -g {input.gtf} {input.quant_dirs} -o {output}'
