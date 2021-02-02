# CDseq

## setting up the workflow

1. Add your raw `.fastq` files to `raw_data/`
2. Update `bin/samples.tsv` to reflect the samples in your dataset.
3. Make sure the information in bin/config.yaml is correct, in particular the genome index, anntaiton file, and genome sequences. 
4. After samples.tsv and config.yaml files are fully setup, load the snakemake module as "module load bbc/snakemake/snakemake-5.28.0". Then, do a dry-run using "snakemake -np". If everything works, the screen will show a list of jobs that the pipeline will run. 
5. When ready to run the pipeline, do "qsub ./bin/run_snake.sh". This will start the job. 

