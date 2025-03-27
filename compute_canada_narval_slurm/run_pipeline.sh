#!/bin/bash
#SBATCH --time=6-23:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --account=def-yeaman

module load nextflow
module load apptainer

source /home/gabnoc/nf-core-env/bin/activate

export NXF_SINGULARITY_CACHEDIR=/project/def-yeaman/NXF_SINGULARITY_CACHEDIR


nextflow run . --ref_genome=/home/gabnoc/scratch/nextflow_stuff/data/ref/GCF_000001735.4_TAIR10.1_genomic.fasta --gff_file=/home/gabnoc/scratch/nextflow_stuff/data/genes/GCF_000001735.4_TAIR10.1_genomic_CHR_MATCHED.gff -profile narval -config nextflow_singularity.config -w /home/gabnoc/scratch/nextflow_stuff/work/
