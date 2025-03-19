#!/usr/bin/env nextflow

params.reads = "./*{1,2}.fastq.gz"
params.outdir = "./"
params.ref_genome = "/lu213/gabriele.nocchi/nextflow_test/GCF_000001735.4_TAIR10.1_genomic.fasta"

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
    tag "Fastp trimming fastq files"
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





process fastaIndex {
    tag "Reference FASTA index building"    
    
    input:
    path reference

    output:
    file "${reference.baseName}.fasta.fai"

    script:
    """
    samtools faidx $reference
    """
}





process gatkIndex {
    tag "Reference GATK index building"

    input:
    path reference

    output:
    file "${reference.baseName}.dict"

    script:
    """
    java -jar /usr/picard/picard.jar CreateSequenceDictionary R=$reference O=${reference.baseName}.dict
    """
}





process bwaIndex {
    tag "Reference BWA index building"    

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
    tag "BWA-mem mapping"
    cpus 2
   
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





process samtoolsSort {
    tag "Samtools sorting bam files"
    cpus 2

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





process addRG {
    tag "Picard add RG to bam files"
    
    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    path "${sorted_bam[0].baseName}_RG.bam"
    script:
    """
    java -jar /usr/picard/picard.jar AddOrReplaceReadGroups -INPUT ${sorted_bam[0]} -OUTPUT ${sorted_bam[0].baseName}_RG.bam -RGID ${sorted_bam[0].baseName} -RGLB ${sorted_bam[0].baseName}_LB -RGPL ILLUMINA -RGPU unit1 -RGSM ${sorted_bam[0].baseName} --VALIDATION_STRINGENCY SILENT
    """
}





process dupRemoval {
    tag "Picard removing duplicates from bam files"

    input:
    path rg_bam

    output:
    path "${rg_bam[0].baseName}_dedup.bam"
    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates -INPUT $rg_bam -OUTPUT ${rg_bam[0].baseName}_dedup.bam -METRICS_FILE ${rg_bam[0].baseName}_DUP_metrics.txt -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT
    """
}





process samtoolsIndex {
    tag "Indexing bam files for indel realignment"

    input:
    path dedup_bam

    output:
    path "${dedup_bam[0].baseName}.bam.bai"
    script:
    """
    samtools index $dedup_bam
    """
}





process realignIndel {
    tag "GATK 3.8 Indel realignment of bam files"

    input:
    path dedup_bam
    path "${dedup_bam[0].baseName}.bam.bai"
    path reference
    file "${reference.baseName}.fasta.fai"
    file "${reference.baseName}.dict"

    output:
    path "${dedup_bam[0].baseName}_realigned.bam"
    script:
    """
    java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $reference -I $dedup_bam -o ${dedup_bam[0].baseName}_intervals.intervals
    java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner -R $reference -I $dedup_bam -targetIntervals ${dedup_bam[0].baseName}_intervals.intervals --consensusDeterminationModel USE_READS  -o ${dedup_bam[0].baseName}_realigned.bam
    """
}




process samtoolsRealignedIndex {
    tag "Indexing bam files for SNP calling"

    input:
    path realigned_bam

    output:
    path "${realigned_bam[0].baseName}.bam.bai"
    script:
    """
    samtools index $realigned_bam
    """
}





process snpCalling {
    tag "SNP calling with bcftools mpileup + call"
    publishDir params.outdir, mode: 'copy'

    input:
    path reference
    file fai
    path bam_files  // This is the collected list of BAM files
    path bai_files

    output:
     path "final_variants_chr_*.vcf.gz"

    script:
    """
     for chr in \$(cut -f1 ${fai}); do
        echo "Processing chromosome: \$chr"
        bcftools mpileup -Ou -f ${reference} -r \$chr ${bam_files} -q 5 -I -a FMT/AD | \
        bcftools call -G - -f GQ -mv -Oz > variants_chr_\${chr}.vcf.gz
        bcftools filter -e 'AC=AN || MQ < 30' variants_chr_\${chr}.vcf.gz -Oz > final_variants_chr_\${chr}.vcf.gz
    done
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
    fai_index = fastaIndex(params.ref_genome)
    gatk_index = gatkIndex(params.ref_genome)
    bwa_index = bwaIndex(params.ref_genome)
    mapped_sam = bwaMap(params.ref_genome, bwa_index, trimmed_ch)
    sorted_bam = samtoolsSort(mapped_sam)
    rg_bam = addRG(sorted_bam)
    dedup_bams = dupRemoval(rg_bam)
    indexed_bams = samtoolsIndex(dedup_bams)
    realigned_bams = realignIndel(dedup_bams,indexed_bams, params.ref_genome, fai_index, gatk_index)
    realigned_bai = samtoolsRealignedIndex(realigned_bams)

    // Extract all BAM paths and collect all of them
    all_bams_ch = realigned_bams.collect()
    all_bai_ch = realigned_bai.collect()
    // Call SNPs
    snpCalling(params.ref_genome,fai_index, all_bams_ch, all_bai_ch)
}

