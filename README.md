# snp_call_nextflow
Learning nextflow by re-implementing a snp calling pipeline



One for docker:

nextflow run snpcall.nf -config nextflow.config --ref_genome /data/home/gabriele/Misc/nextflow/nextflow_docker/GCF_000001735.4_TAIR10.1_genomic.fasta

One for apptainer/singularity:

nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /data/home/gabriele/Misc/nextflow/nextflow_singularity/GCF_000001735.4_TAIR10.1_genomic.fasta --gff_file /data/home/gabriele/Misc/nextflow/nextflow_singularity/genes.gff



ALWAYS GIVE FULL PATH FOR REFERENCE AND GFF FILE

REFERENCE MUST HAVE .fasta SUFFIX (change it to .fasta if .fa)

GFF FILE MUST HAVE .gff SUFFIX
