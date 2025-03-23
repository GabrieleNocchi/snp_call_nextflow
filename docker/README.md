This pipeline takes fastq reads and a reference genome and will produce a final minimally filtered vcf (removing SNPs where all indidivuals are homozyogous ALT and any SNP with MQ < 30).

This pipeline works the same ways as the singularity/apptainer based pipeline, however it does not output depth of coverage statistics per sample. 

Therefore, it does not have the gff_file parameter (--gff_file)    

CD to the work dir where your fastq are and run:  

<pre>nextflow run snpcall_singularity.nf -config nextflow_singularity.config --ref_genome /path/to/reference_genome.fasta</pre>


<b>Available options:</b>

--reads (default CWD: ./*{1,2}.fastq.gz)

--outdir (default CWD: ./)

--ref_genome (No default, give full path)


<b>IMPORTANT</b>

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix

Need to pull all the docker images and link them to the directory where you save them in the config file nextflow.config
