# snp_call_nextflow
Learning nextflow by re-implementing a snp calling pipeline. 

This pipeline takes fastq reads, a reference genome and a genes gff file and will produce a final minimally filtered vcf (removing SNPs where all indidivuals are homozyogous ALT and any SNP with MQ < 30) and 3 depth statistics files per dataset (genes depth , windows depth, wg depth)



To run the apptainer/singularity based pipeline CD to the work dir where your fastq are and run:

nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /path/to/reference_genome.fasta --gff_file /path/to/genes.gff


Available options:

--reads (default CWD: ./*{1,2}.fastq.gz)

--outdir (default CWD: ./)

--ref_genome (No default, give full path)

--gff_file (No default, give full path)


IMPORTANT

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix
