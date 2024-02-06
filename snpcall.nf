#!/usr/bin/env nextflow

params.reads = "/lu213/gabriele.nocchi/*_{1,2}.fastq.gz"
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
    tag "Fastp trimming"
    cpus 2
  
    
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
    tag "BWA index building"    
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
    tag "BWA mapping"
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
    path "${sample_id}.sam"

    script:
    """
    bwa mem -t $task.cpus $reference ${trimmed_reads[0]} ${trimmed_reads[1]} > ${sample_id}.sam
    """
}


process samSort {
    tag "Samtools sorting bams"
    cpus 4

    input:
    path sample_sam

    output:
    tuple val("${sample_sam.baseName}"), path("${sample_sam.baseName}{_sorted.bam,_sorted.bam.bai}")
    
    script:
    """
    samtools view -Sb -q 10 $sample_sam > temp666.bam
    rm $sample_sam
    samtools sort -n -o temp777.bam temp666.bam
    samtools fixmate -m temp777.bam temp888.bam
    samtools sort --threads $task.cpus temp888.bam > ${sample_sam.baseName}_sorted.bam
    rm temp666.bam temp777.bam temp888.bam
    samtools index ${sample_sam.baseName}_sorted.bam
    """
}



process dupRemoval {
    tag "Samtools removing duplicates"
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val("${sorted_bam[0].baseName}"),path("${sorted_bam[0].baseName}{_dedup.bam,_dedup.bam.bai}")

    script:
    """
    samtools markdup -r -s ${sorted_bam[0]} ${sorted_bam[0].baseName}_dedup.bam
    samtools index ${sorted_bam[0].baseName}_dedup.bam
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
    mapped_sam = bwaMap(params.ref_genome, bwa_index, trimmed_ch)
    sorted_bam = samSort(mapped_sam)
    dupRemoval(sorted_bam)
}

