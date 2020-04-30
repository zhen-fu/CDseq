#PBS -l walltime=100:00:00
#PBS -l mem=128gb
#PBS -m ae
#PBS -N trim_2

cd /secondary/projects/bbc/research/PFEG_20200427_CDseq
module load bbc/trim_galore/trim_galore-0.6.0
module load bbc/pigz/pigz-2.4

trim_galore \
--paired \
/secondary/projects/bbc/research/PFEG_20200427_CDseq/raw_data/SE6114_SA58949_S2_L005_R1_001.fastq.gz /secondary/projects/bbc/research/PFEG_20200427_CDseq/raw_data/SE6114_SA58949_S2_L005_R2_001.fastq.gz \
--output_dir analysis/trim_test \
--cores 4 \
-q 20 \
--fastqc \
--hardtrim5 50
