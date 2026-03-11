This is a pipeline to process scNMT sequencing data, intended for DKFZ A290 internal use.



# Important notes

- A lot of steps of the pipeline won't work on a Windows volume (e.g. the O drive). As a result, **the pipeline folder, the output directory and the reference files all need to reside on a Linux volume** (e.g. scratch). In addition, it is recommended to have no white space in any folder name or file name.
- It is highly recommeded to use absolute paths when setting up the config file (see below).
- STAR requires around 30 GB of RAM already for the human or mouse genome (~3 GB). Make sure your machine has enough RAM left.

# Environment setup

Running the pipeline requires having the following commands available from command line.

- `conda`
- `snakemake`

One way to get `conda` is to install Miniforge, following instructions on [their GitHub repo](https://github.com/conda-forge/miniforge), but feel free to choose any alternative.

Once `conda` is available, create a conda environment called "snakemake" with `conda create -n snakemake`. This will be our main environment.

Now you can activate the environment and install the required packages:

```bash
conda activate snakemake
conda install snakemake pandas<2.0 xlrd openpyxl
```

Alternatively, you can use the the provided YAML file

```bash
conda env create -n snakemake -f envs/snakemake.yaml
```

> `pandas>=2.0` has compatibility issues so we need a lower version. `xlrd` and `openpyxl` are needed only when your metadata are `.xls` or `.xlsx`, respectively. If that is not the case, feel free to drop these.

# Reference genome and annotations

Both the reference genome and gene annotations can be downloaded from [EMSELBL](https://www.ensembl.org/index.html).

For mouse GRCm39 of release 115, use the following commands:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz
```

Similarly, for human GRCh38 release 115, use:

```bash
curl -O https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
curl -O https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz
```

These reference files are typically gzipped (i.e. ending with `.gz`) and must be unzipped before running the pipeline. This can be done with e.g. `gzip -d [filename]` (or `gzip -dk [filename]` if you want to keep the unzipped file).

# Configuration

The config file `snakeconfig.yaml` is where you configure the pipeline parameters.

```yaml
dataset: [dataset_name]
outdir: [path_to_output_directory]

ref_genome:
  seq: [path_to_reference_genome_fasta_file]
  genes: [path_to_annotation_gtf_file]

ilse_info:
  [ILSe ID]:
    metadata: [path_to_metadata_file]
    fastqdir: [path_to_fastq_directory]
```

- `[dataset_name]`: a name given to the dataset. This will be used to name the output files, but it shouldn't have any effect on the pipeline itself.
- `[path_to_output_directory]`: the directory where all outputs and log files will go.
- `[path_to_reference_genome_fasta_file]`, `[path_to_annotation_gtf_file]`: path to the reference genome and gene annotation files. Again, these must be unzipped (see section above).
- `[path_to_metadata_file]`: metadata file downloaded from ILSe website, usually in the format `[ILSe ID]-result.xls`. Remember to remove rows that do not belong to the dataset of interest.
- `[path_to_fastq_directory]`: directory to find the raw sequenceing results, usually on the NGS drive.

> You can add any number of ILSe IDs to the `ilse_info` section, each with a metadata field and a fastq field. The pipeline will identify samples by the "Sample Name" column in the metadata file, so if the same sample name appears in two different metadata files, the pipeline will treat them as one sample being sequenced twice. Therefore, make sure that sample names are consistent across sequencing runs.

It is highly recommended to use absolute path for all entries in the config file, since relative paths can lead to unexpected errors.


# Running the pipeline

Change to the directory containing the script and the activate the snakemake environment.

```bash
cd /path/to/pipeline/directory
conda activate snakemake
```

Use the following command to run the pipeline.

```bash
snakemake --cores 24 --sdm conda
```

Change the number of cores to match your system's capabilities. For a dry run, add the `-n` flag to the end of the command.


# Known issues

If STAR is interrupted from outside (e.g. command line interruption or killed due to memory shortage), the loaded genome might stay in memory and become locked, preventing the next STAR instance to load any genome at all.
When this happens, use `STAR --genomeDir [/path/to/STAR/indices] --genomeLoad Remove`.
