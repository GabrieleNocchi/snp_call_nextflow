#!/usr/bin/env nextflow

params.reads = "/lu213/gabriele.nocchi/*_R{1,2}*.gz"
params.outdir = "./"
params.ref_genome = "/lu213/gabriele.nocchi/GCF_000001735.4_TAIR10.1_genomic.fasta"

log.info """\
    S N P   C A L L I N G   P I P E L I N E
    ===================================
    reference    : ${params.ref_genome}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)


/*
 * Processes
 */
process trimSequences {
    tag "Fastp on $sample_id"
    cpus 4
  
    
    input:
    tuple val(sample_id), path(reads)

    
    output:
    tuple val(sample_id), path("${sample_id}_{1,2}_trimmed.fastq.gz")

    script:
    """  
    fastp -w $task.cpus -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_1_trimmed.fastq.gz -O ${sample_id}_2_trimmed.fastq.gz
    """

}



process bwaIndex {
        
    input:
    path reference
    
    output:
    file "${reference.baseName}.fasta.amb"
    file "${reference.baseName}.fasta.ann"
    file "${reference.baseName}.fasta.bwt"
    file "${reference.baseName}.fasta.pac"
    file "${reference.baseName}.fasta.sa"    

    script:
    """
    bwa index -a bwtsw $reference 
    """
}

process bwaMap {
   
    publishDir params.outdir, mode: 'copy'
    cpus 4
   
    input:
    path reference
    file "${reference.baseName}.fasta.amb"
    file "${reference.baseName}.fasta.ann"
    file "${reference.baseName}.fasta.bwt"
    file "${reference.baseName}.fasta.pac"
    file "${reference.baseName}.fasta.sa"
    tuple val(sample_id), path(trimmed_reads)

    output:
    file "final.sam"

    script:
    """
    bwa mem -t $task.cpus $reference ${trimmed_reads[0]} ${trimmed_reads[1]} > final.sam
    """
}



/*
 * Define the workflow
 */
workflow {

    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
   
   
    trimmed_ch = trimSequences(read_pairs_ch)
    bwa_index = bwaIndex(params.ref_genome)
    bwaMap(params.ref_genome, bwa_index, trimmed_ch)
}

