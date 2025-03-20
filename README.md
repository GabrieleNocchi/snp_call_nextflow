# snp_call_nextflow
Learning nextflow by re-implementing a snp calling pipeline. 

This pipeline takes fastq reads, a reference genome and a genes gff file and will produce a final minimally filtered vcf (removing SNPs where all indidivuals are homozyogous ALT and any SNP with MQ < 30) and 3 depth statistics files per dataset (genes depth , windows depth, wg depth)



Apptainer/singularity based pipeline (default your fastq must be in the work dir, output will go in the work dir):

nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /data/home/gabriele/Misc/nextflow/nextflow_singularity/GCF_000001735.4_TAIR10.1_genomic.fasta --gff_file /data/home/gabriele/Misc/nextflow/nextflow_singularity/genes.gff


Available options:

--reads (default: ./*{1,2}.fastq.gz)

--outdir (default: ./)

--ref_genome (No default, give full path)

--gff_file (No default, give full path)


IMPORTANT

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix
