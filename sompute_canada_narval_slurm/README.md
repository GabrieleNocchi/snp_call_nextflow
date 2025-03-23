Run the singularity/apptainer based pipeline in compute canada narval



From a login node:

<pre>module purge # Make sure that previously loaded modules are not polluting the installation 
module load python/3.11
module load rust # New nf-core installations will err out if rust hasn't been loaded
module load postgresql # Will not use PostgresSQL here, but some Python modules which list psycopg2 as a dependency in the installation would crash without it.
python -m venv nf-core-env
source nf-core-env/bin/activate
python -m pip install nf_core==2.13</pre>


Now, create or edit the file (you probably have to create it):  ~/.nextflow/config
You can use the config below, but change def-yeaman to your account name

<pre>params {
    config_profile_description = 'Alliance HPC config'
    config_profile_contact = 'https://docs.alliancecan.ca/wiki/Technical_support'
    config_profile_url = 'docs.alliancecan.ca/wiki/Nextflow'
}


singularity {
  enabled = true
  autoMounts = true
}

apptainer {
  autoMounts = true
}

process {
  executor = 'slurm'
  clusterOptions = '--account=def-yeaman'
  maxRetries = 1
  errorStrategy = { task.exitStatus in [125,139] ? 'retry' : 'finish' }
  memory = '4GB'
  cpu = 1
  time = '3h'
}

executor {
  pollInterval = '60 sec'
  submitRateLimit = '60/1min'
  queueSize = 100
}

profiles {
  beluga {
    max_memory='186G'
    max_cpu=40
    max_time='168h'
  }
  narval {
    max_memory='249G'
    max_cpu=64
    max_time='168h'
  }
}
</pre>



Now, download all the singularity images needed: https://github.com/RepAdapt/singularity/blob/main/RepAdaptSingularity.md

Place them here:
<pre>mkdir /project/<def-group>/NXF_SINGULARITY_CACHEDIR</pre>
then:
<pre>export NXF_SINGULARITY_CACHEDIR=/project/<def-group>/NXF_SINGULARITY_CACHEDIR</pre>

Also add the above export command to your ~/.bashrc


Now we have everything ready to start the workflow.
To do so, start a screen session:

<pre>screen -S snp_calling
module load nextflow
module load apptainer

source ~/nf-core-env/bin/activate
nextflow run . --ref_genome=/full/path/to/data/ref/GCF_000001735.4_TAIR10.1_genomic.fasta --gff_file=/full/path/to/data/genes/genes.gff -profile narval -config nextflow_singularity.config -w /home/gabnoc/scratch/nextflow_stuff/work/
</pre>


<b>IMPORTANT</b>

--reads (default CWD: "/home/gabnoc/scratch/nextflow_stuff/data/*{1,2}.fastq.gz") ### this can be changed by passing --reads as flag to the command, to match your reads location path and names patterns (ie. fq.gz)

--outdir (default CWD: "/home/gabnoc/scratch/nextflow_stuff/output/") ### this can be changed to any directory

--ref_genome (No default, give full path)

--gff_file (No default, give full path)


<b>IMPORTANT</b>

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix

Need to pull all the apptainer images as sif files and link them to the directory where you save them in the config file nextflow_singularity.config
