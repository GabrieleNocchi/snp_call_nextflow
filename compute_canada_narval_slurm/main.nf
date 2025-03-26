#!/usr/bin/env nextflow

params.reads = "/home/gabnoc/scratch/nextflow_stuff/data/*{1,2}.fastq.gz"
params.outdir = "/home/gabnoc/scratch/nextflow_stuff/output/"
params.ref_genome = ""
params.gff_file = ""


log.info """\
    S N P   C A L L I N G   P I P E L I N E
    ===================================
    reference    : ${params.ref_genome}
    reads        : ${params.reads}
    gff          : ${params.gff_file}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)



/*
 * Processes
 */
process trimSequences {
    tag "Fastp trimming fastq files"
    cpus 4
    memory '4GB'
    time '6h'
    errorStrategy 'ignore'
  
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
    cpus 1
    memory '4GB'
    time '6h'

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
    cpus 1
    memory '4GB'
    time '6h'

    input:
    path reference

    output:
    file "${reference.baseName}.dict"

    script:
    """
    picard CreateSequenceDictionary R=$reference O=${reference.baseName}.dict
    """
}





process bwaIndex {
    tag "Reference BWA index building"    
    cpus 1
    memory '4GB'
    time '6h'

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
    cpus 4
    memory '4GB'
    time '6h'
    errorStrategy 'ignore'
   
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
    cpus 4
    memory '32GB'
    time '12h'
    errorStrategy 'ignore'

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
    cpus 1
    memory '4GB'
    time '6h'
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    path "${sorted_bam[0].baseName}_RG.bam"
    script:
    """
    picard AddOrReplaceReadGroups -INPUT ${sorted_bam[0]} -OUTPUT ${sorted_bam[0].baseName}_RG.bam -RGID ${sorted_bam[0].baseName} -RGLB ${sorted_bam[0].baseName}_LB -RGPL ILLUMINA -RGPU unit1 -RGSM ${sorted_bam[0].baseName} --VALIDATION_STRINGENCY SILENT
    """
}





process dupRemoval {
    tag "Picard removing duplicates from bam files"
    cpus 1
    memory '4GB'
    time '6h'
    errorStrategy 'ignore'

    input:
    path rg_bam

    output:
    path "${rg_bam[0].baseName}_dedup.bam"
    script:
    """
    picard MarkDuplicates -INPUT $rg_bam -OUTPUT ${rg_bam[0].baseName}_dedup.bam -METRICS_FILE ${rg_bam[0].baseName}_DUP_metrics.txt -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT
    """
}






process realignIndel {
    tag "GATK 3.8 Indel realignment of bam files"
    cpus 1
    memory '64GB'
    time '48h'
    errorStrategy 'ignore'

    input:
    path dedup_bam
    path reference
    file "${reference.baseName}.fasta.fai"
    file "${reference.baseName}.dict"

    output:
    path "${dedup_bam[0].baseName}_realigned.bam"
    script:
    """
    gatk3 -T RealignerTargetCreator -R $reference -I $dedup_bam -o ${dedup_bam[0].baseName}_intervals.intervals -U ALLOW_UNINDEXED_BAM
    gatk3 -T IndelRealigner -R $reference -I $dedup_bam -targetIntervals ${dedup_bam[0].baseName}_intervals.intervals --consensusDeterminationModel USE_READS  -o ${dedup_bam[0].baseName}_realigned.bam -U ALLOW_UNINDEXED_BAM
    """
}




process samtoolsRealignedIndex {
    tag "Indexing bam files for SNP calling"
    cpus 1
    memory '4GB'
    time '6h'
    errorStrategy 'ignore'

    input:
    path realigned_bam

    output:
    path "${realigned_bam[0].baseName}.bam.bai"
    script:
    """
    samtools index $realigned_bam
    """
}




process prepareDepth {
    tag "Preparing files for depth statistics"
    cpus 1
    memory '4GB'
    time '6h'
    publishDir params.outdir, mode: 'copy'

    input:
    path reference_fai
    path gff_file

    output:
    file "genome.bed"
    file "windows.bed"
    file "windows.list"
    file "genes.bed"
    file "genes.list"
   
    script:
    """
    awk '{print \$1"\\t"\$2}' $reference_fai > genome.bed

    awk -v w=5000 '{chr = \$1; chr_len = \$2;
    for (start = 0; start < chr_len; start += w) {
        end = ((start + w) < chr_len ? (start + w) : chr_len);
        print chr "\\t" start "\\t" end;
       }
    }' $reference_fai > windows.bed

    awk -F "\\t" '{print \$1":"\$2"-"\$3}' windows.bed | sort -k1,1 > windows.list
    awk '\$3 == "gene" {print \$1"\\t"\$4"\\t"\$5}' $gff_file | uniq > genes.bed

    cut -f1 $reference_fai | while read chr; do awk -v chr=\$chr '\$1 == chr {print \$0}' genes.bed | sort -k2,2n; done > genes.sorted.bed
    mv genes.sorted.bed genes.bed

    awk -F "\\t" '{print \$1":"\$2"-"\$3}' genes.bed | sort -k1,1 > genes.list
    """
}




process calculateDepth {
    tag "Extracting depth statistics from bam files"
    cpus 1
    memory '16GB'
    time '6h'

    input:
    path realigned_bam

    output:
    path "${realigned_bam[0].baseName}.depth"
    script:
    """
    samtools depth -aa $realigned_bam > ${realigned_bam[0].baseName}.depth
    """
}




process calculateGenesDepth {
    tag "Extracting genes depth statistics from bam files"
    cpus 1
    memory '32GB'
    time '6h'

    input:
    path depth_file
    file "genome.bed"
    file "windows.bed"
    file "windows.list"
    file "genes.bed"
    file "genes.list"

    output:
    path "${depth_file[0].baseName}-genes.sorted.tsv"
    

    script:
    """
    awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3}' $depth_file | bedtools map -a genes.bed -b stdin -c 4 -o mean -null 0 -g genome.bed | awk -F "\\t" '{print \$1":"\$2"-"\$3"\\t"\$4}' | sort -k1,1 > ${depth_file[0].baseName}-genes.tsv
    awk 'NR==FNR{a[\$1]=\$0; next} \$1 in a{print a[\$1]"\\t"\$2}' genes.list ${depth_file[0].baseName}-genes.tsv > ${depth_file[0].baseName}-genes.sorted.tsv
    """
}




process calculateWindowsDepth {
    tag "Extracting windows depth statistics from bam files"
    cpus 1
    memory '32GB'
    time '6h'

    input:
    path depth_file
    file "genome.bed"
    file "windows.bed"
    file "windows.list"
    file "genes.bed"
    file "genes.list"

    output:
    path "${depth_file[0].baseName}-windows.sorted.tsv"


    script:
    """
    awk '{print \$1"\\t"\$2"\\t"\$2"\\t"\$3}' $depth_file | bedtools map -a windows.bed -b stdin -c 4 -o mean -null 0 -g genome.bed | awk -F "\\t" '{print \$1":"\$2"-"\$3"\\t"\$4}' | sort -k1,1 > ${depth_file[0].baseName}-windows.tsv
    awk 'NR==FNR{a[\$1]=\$0; next} \$1 in a{print a[\$1]"\\t"\$2}' windows.list ${depth_file[0].baseName}-windows.tsv > ${depth_file[0].baseName}-windows.sorted.tsv
    """
}




process calculateWgDepth {
    tag "Extracting WG depth statistics from bam files"
    cpus 1
    memory '32GB'
    time '6h'

    input:
    path depth_file
    file "genome.bed"
    file "windows.bed"
    file "windows.list"
    file "genes.bed"
    file "genes.list"

    output:
    path "${depth_file[0].baseName}-wg.txt"


    script:
    """
    awk '{sum += \$3; count++} END {if (count > 0) print sum/count; else print "No data"}' $depth_file > ${depth_file[0].baseName}-wg.txt
    """
}





process joinDepth {
    tag "Joining depth statistics across samples"
    cpus 1
    memory '8GB'
    time '6h'
    publishDir params.outdir, mode: 'copy'

    input:
    path bam_files
    path gene_statistics_files
    path window_statistics_files
    path wg_statistics_files    

    output:
    path "combined_windows.tsv"
    path "combined_genes.tsv"
    path "combined_wg.tsv"
    
    script:
    """
    echo -e "${bam_files.collect { it.baseName.replace('.bam', '') }.join('\\n')}" > list.txt
    echo -e "location\\t\$(cut -f2 list.txt | sort | uniq | paste -s -d '\\t')" > depthheader.txt
    cut -f2 list.txt | sort | uniq > samples.txt

    while read samp; do cut -f2 \${samp}-windows.sorted.tsv > \${samp}-windows.sorted.depthcol ; done < samples.txt
    paste \$(sed 's/^/.\\//' samples.txt | sed 's/\$/-windows.sorted.tsv/' | head -n 1) \$(sed 's/^/.\\//' samples.txt | sed 's/\$/-windows.sorted.depthcol/' | tail -n +2) > combined-windows.temp
    cat depthheader.txt combined-windows.temp > combined_windows.tsv

    while read samp; do cut -f2 \${samp}-genes.sorted.tsv > \${samp}-genes.sorted.depthcol ; done < samples.txt
    paste \$(sed 's/^/.\\//' samples.txt | sed 's/\$/-genes.sorted.tsv/' | head -n 1) \$(sed 's/^/.\\//' samples.txt | sed 's/\$/-genes.sorted.depthcol/' | tail -n +2) > combined-genes.temp
    cat depthheader.txt combined-genes.temp > combined_genes.tsv

    while read samp; do echo -e \$samp"\\t"\$(cat ./\$samp-wg.txt); done < samples.txt > combined_wg.tsv
    """
}





process snpCalling {
    tag "SNP calling with bcftools mpileup + call"
    cpus 1
    memory '16GB'
    time '48h'

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




process concatVCFs {
    tag "Concatenating chromosomes VCFs"
    cpus 1
    memory '16GB'
    time '23h'
    publishDir params.outdir, mode: 'copy'

    input:
    path vcfs

    output:
    path "final_variants.vcf.gz"
    path "final_variants.vcf.gz.tbi"

    script:
    """
    bcftools concat $vcfs -Oz > final_variants.vcf.gz
    tabix -p vcf final_variants.vcf.gz
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
    depth_prep_files = prepareDepth(fai_index, params.gff_file)
    gatk_index = gatkIndex(params.ref_genome)
    bwa_index = bwaIndex(params.ref_genome)
    mapped_sam = bwaMap(params.ref_genome, bwa_index, trimmed_ch)
    sorted_bam = samtoolsSort(mapped_sam)
    rg_bam = addRG(sorted_bam)
    dedup_bams = dupRemoval(rg_bam)
    realigned_bams = realignIndel(dedup_bams, params.ref_genome, fai_index, gatk_index)
    realigned_bai = samtoolsRealignedIndex(realigned_bams)
    depth_file = calculateDepth(realigned_bams)
    depth_stats_genes = calculateGenesDepth(depth_file, depth_prep_files)
    depth_stats_windows = calculateWindowsDepth(depth_file, depth_prep_files)
    depth_stats_wg = calculateWgDepth(depth_file, depth_prep_files)

    // Extract all BAM paths and collect all of them
    all_bams_ch = realigned_bams.collect()
    all_bai_ch = realigned_bai.collect()


    // Combine and map depth stats output
    all_genes_ch = depth_stats_genes.collect()
    all_windows_ch = depth_stats_windows.collect()
    all_wg_ch = depth_stats_wg.collect()
    joinDepth(all_bams_ch, all_genes_ch, all_windows_ch, all_wg_ch)
    

    // Call SNPs
   vcfs = snpCalling(params.ref_genome,fai_index, all_bams_ch, all_bai_ch)
   all_vcfs_ch = vcfs.collect()
   concatVCFs(all_vcfs_ch)
}

