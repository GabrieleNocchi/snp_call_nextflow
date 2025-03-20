# snp_call_nextflow
Learning nextflow by re-implementing a snp calling pipeline. 

This pipeline takes fastq reads, a reference genome and a genes gff file and will produce a final minimally filtered vcf (removing SNPs where all indidivuals are homozyogous ALT and any SNP with MQ < 30) and 3 depth statistics files per dataset (genes depth , windows depth, wg depth)


To run in the directory where your raw paired end fastq files are



Apptainer/singularity based pipeline:

nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /data/home/gabriele/Misc/nextflow/nextflow_singularity/GCF_000001735.4_TAIR10.1_genomic.fasta --gff_file /data/home/gabriele/Misc/nextflow/nextflow_singularity/genes.gff



IMPORTANT

ALWAYS GIVE FULL PATH FOR REFERENCE AND GFF FILE 

REFERENCE MUST HAVE .fasta SUFFIX (change it to .fasta if yours is .fa)

GFF FILE MUST HAVE .gff SUFFIX

PAIRED FASTQ FILES MUST BE IN WORK DIR AND FOLLOW THIS NAMING: *{1,2}.fastq.gz 

else this regex can be changed and given as a flag --reads to match your naming
