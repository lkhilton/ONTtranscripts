# ONTtranscripts

## Requirements

Conda: `https://docs.conda.io/en/latest/miniconda.html`

Snakemake: `conda install -c bioconda snakemake`

## Data

Download and unzip the tarballs containing fast5 and fastq files to the data directory
```
data
├── fast5
│   ├── IX7311
│   └── IX7312
├── fastq
└── metadata
```

## Reference files

Place or symlink a whole genome fasta file in `ref`. I used an hg19 build. 

Generate a minimap index file: `minimap2 -d GRCh37-lite.mmi GRCh37-lite.fa`. 

The transcript reference fasta and fai index files are included in this repo. 

## Config

Modify `config/config.yaml` to include local paths to your reference files. 


## Running Snakemake

Once the data directory is populated and the reference and config files are updated, `cd` into the top level directory and launch with: 
```
snakemake -j {num_cpus} --use-conda --latency-wait 120
```
`{num_cpus}`: INT > 8
