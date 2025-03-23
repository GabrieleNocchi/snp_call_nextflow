This pipeline takes fastq reads, a reference genome and a genes gff file and will produce a final minimally filtered vcf (removing SNPs where all indidivuals are homozyogous ALT and any SNP with MQ < 30) and 3 depth statistics files per dataset (genes depth , windows depth, wg depth)

To run the singularity/apptainer based pipeline CD to the work dir where your fastq are and run:
<pre>nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /path/to/reference_genome.fasta --gff_file /path/to/genes.gff</pre>


<b>Available options:</b>

--reads (default CWD: ./*{1,2}.fastq.gz) ### this can be changed by passing --reads as flag to the command, to match your reads location path and names patterns

--outdir (default CWD: ./) ### this can be changed to any directory

--ref_genome (No default, give full path)

--gff_file (No default, give full path)


<b>IMPORTANT</b>

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix

Need to pull all the apptainer images as sif files and link them to the directory where you save them in the config file nextflow_singularity.config
