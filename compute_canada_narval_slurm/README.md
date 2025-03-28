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

This is like a general config for EVERY workflow you will run with nextflow in this cluster (narval). You can copy and paste the text below into it, but change def-group to your account name



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
  clusterOptions = '--account=def-group'
  maxRetries = 1
  errorStrategy = { task.exitStatus in [125,139] ? 'retry' : 'finish' }
  memory = '4GB'
  cpu = 1
  time = '1h'
}

executor {
  pollInterval = '60 sec'
  submitRateLimit = '60/1min'
  queueSize = 900
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

Do not worry about the cpus, memory and time of the slurm process. These slurm global options will be overwritten by the cpus, memory and time specified for each process of the workflow in main.nf


Now, download all the singularity images needed: https://github.com/RepAdapt/singularity/blob/main/RepAdaptSingularity.md

Place them here (replace def-group with your account name):
<pre>mkdir /project/def-group/NXF_SINGULARITY_CACHEDIR</pre>
then:
<pre>export NXF_SINGULARITY_CACHEDIR=/project/def-group/NXF_SINGULARITY_CACHEDIR</pre>

Also add the above export command to your ~/.bashrc


Now we have everything ready to start the workflow.

The workflow can be run on a computing node, using the script <b>run_pipeline.sh</b> (submit with sbatch -- give this job only 1 cpu, 4 GB RAM but MAX available run time which is 7 days in narval and beluga). 

<b>AVAILABLE OPTIONS</b>

--reads (default CWD: "/home/gabnoc/scratch/nextflow_stuff/data/*{1,2}.fastq.gz") ### this can be changed by passing --reads as flag to the command, to match your reads location path and names patterns (ie. fq.gz)

--outdir (default CWD: "/home/gabnoc/scratch/nextflow_stuff/output/") ### this can be changed to any directory with flag --outdir

--ref_genome (No default, give full path)

--gff_file (No default, give full path)


<b>IMPORTANT</b>

The reference genome file MUST HAVE .fasta suffix (change it to .fasta if yours is .fa)

The GFF file MUST HAVE .gff suffix

Need to pull all the apptainer images as sif files and link them to the directory where you saved  them (/project/def-group/NXF_SINGULARITY_CACHEDIR) in the config file nextflow_singularity.config
