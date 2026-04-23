# version control
from snakemake.utils import min_version
min_version('9.0')

# Python libraries to be used
from collections import defaultdict
from os.path import dirname, basename, join
import re

import pandas as pd
import yaml


# read config and overwrite default
configfile: 'snakeconfig.yaml'
with open('defaultconfig.yaml') as f:
    default_config = yaml.safe_load(f)
config = default_config | config

# sanity checks
for section in ('dataset', 'reference', 'outdir', 'ilse_info'):
    if section not in config:
        raise ValueError(f'No "{section}" section found in config file')

if config['pipeline'] not in ('star_umite', 'salmon', 'biscuit_methscan'):
    raise ValueError(f'Unexpected pipeline {config['pipeline']}')

# for subsection in ('genome', 'genes'):
#     if not subsection in config['reference']:
#         raise ValueError(f'No "reference - {subsection}" found in config file')

# if config['pipeline'] == 'salmon' and 'transcriptome' not in config['reference']:
#     raise ValueError('Salmon pipeline requires "reference - transcriptome" in config file')

config['ilse_info']['metadata'] = config['ilse_info']['metadata'].split()
config['ilse_info']['fastqdir'] = config['ilse_info']['fastqdir'].split()

if not config['ilse_info']['metadata']:
    raise ValueError('No metadat file path provided')
if not config['ilse_info']['metadata']:
    raise ValueError('No FASTQ search directory provided')

paths_to_check = config['ilse_info']['metadata'] + config['ilse_info']['fastqdir']
for path in paths_to_check:
    if not os.path.isabs(path):
        raise ValueError('Must use absolute paths in the config file')
    elif not os.path.exists(path):
        raise ValueError(f'Path {path} does not exist')

# construct objects needed for pipeline
sample_to_fqid = defaultdict(list)
fqid_to_dir = {}

for metadata in config['ilse_info']['metadata']:
    metadata_ext = os.path.splitext(metadata)[-1]
    if metadata_ext == '.xls' or metadata_ext == '.xlsx':
        df_ilse = pd.read_excel(metadata, index_col='Sample Name', dtype=str)
    elif metadata_ext == '.csv':
        df_ilse = pd.read_csv(metadata, index_col='Sample Name', dtype=str)
    else:
        raise ValueError(f'Unexpected metadata extension "{metadata_ext}".')

    for sample, row in df_ilse.iterrows():
        fqid = row['Unique ID / Lane']
        found = False
        for searchdir in config['ilse_info']['fastqdir']:
            if fqid in os.listdir(searchdir):
                fqid_to_dir[fqid] = searchdir
                # sample_to_read1[sample].append(join(searchdir, f'{fqid}/fastq/{fqid}_R1.fastq.gz'))
                # sample_to_read2[sample].append(join(searchdir, f'{fqid}/fastq/{fqid}_R2.fastq.gz'))
                sample_to_fqid[sample].append(fqid)
                found = True
                break # If we have repeated FASTQ ID, only the first one will be used
        if not found:
            raise ValueError(f'Could not find FASTQ ID {fqid} for sample {sample} in any of the specified directories') 




# fetch rules according to specified pipeline
include: f'modules/{config['pipeline']}.smk'
# if config['pipeline'] == 'star_umite':
#     include: 'modules/star_umite.smk'
# elif config['pipeline'] == 'salmon':
#     include: 'modules/salmon.smk'
# elif config['pipeline'] == 'biscuit_methscan':
#     include: 'modules/biscuit_methscan.smk'
# else:
#     raise ValueError(f'Invalid pipeline name "{config['pipeline']}"')


# constrain wildcards
wildcard_constraints:
    sample = r'\w+'


# define output
rule all:
    default_target: True
    input:
        #expand(rules.extract_methylation.output, sample='A1')
        join(config['outdir'], f'{config['dataset']}.{config['pipeline']}.h5ad')
    run:
        # if pipeline successfully finishes, write config to outdir
        from datetime import datetime
        with open(join(config['outdir'], 'last_success_config.yaml'), 'w') as f:
            f.write(
                f'# These are configurations used for the most recent successful run of the pipeline, which was {datetime.today().isoformat(sep=" ", timespec="seconds")}.\n'
            )
            yaml.dump(config, f, sort_keys=False)