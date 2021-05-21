#!/usr/bin/env nextflow

// Used CLIP-seq nfcore pipeline as a template, which scripts to run, etc
// https://nf-co.re/clipseq/
// Assuming inputs are gzipped fasta and gtf files, here hardcoded in


// Inputs, Test Dataset plus relevant ref genome
// Oisin


params.outdir = "./"
fasta_path = ("/data/MSc/2021/clipseq/sirna_trimmed_chr20.fa")
fastq_path = ("/data/MSc/2021/clipseq/sirna_trimmed_chr20.fq.gz")
gtf_path = ("/data/MSc/2021/clipseq/ref/chr20.gtf")
refgenome_path = ("/data/MSc/2021/clipseq/ref/chr20.fa")
smrna_path = ("/data/MSc/2021/clipseq/Homo_sapiens.smallRNA.fa")
star_path = ("/data/MSc/2021/clipseq/index_star_ch20")
multiqc_config = ("/data/omccaffrey/MA5112/nextflow_assignment/running_clipseq.nf/config_and_docs/multiqc_config.yaml")
output_docs = ("/data/omccaffrey/MA5112/nextflow_assignment/running_clipseq.nf/config_and_docs/output.md")

// Setup channels
ch_smrna_fasta = Channel.value(smrna_path)
ch_ref_fai = Channel.value(refgenome_path)
ch_ref = Channel.value(refgenome_path)
ch_gtf_star = Channel.value(gtf_path)
ch_fastq_fastqc_pretrim = Channel.value(fastq_path)
ch_star_index = Channel.value(star_path)
ch_fasta = Channel.value(refgenome_path)
ch_fastq_umi = Channel.value(fastq_path)
ch_fasta_pureclip = Channel.value(refgenome_path)
ch_multiqc_config = Channel.value(multiqc_config)
ch_output_docs = Channel.value(output_docs)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()



/*
 * Generating premapping index
 */

process generate_premap_index{

	tag "$smrna_fasta"

    	input:
    	path(smrna_fasta) from ch_smrna_fasta

    	output:

    	path("${smrna_fasta.simpleName}.*.bt2") into ch_bt2_index

    	script:
    	"""
    	bowtie2-build --threads $task.cpus $smrna_fasta ${smrna_fasta.simpleName}
    	"""
}

// Generate fai from reference file
process generate_fai{

  tag "$ref"

    	input:
   		path(ref) from ch_ref_fai

   		output:
   		path("*.fai") into (ch_fai_crosslinks, ch_fai_icount, ch_fai_icount_motif, ch_fai_paraclu_motif, ch_fai_pureclip_motif, ch_fai_piranha_motif)

    	script:
    	"""
    	samtools faidx $ref
    	"""
}

// Generate STAR index

process star_index{

				script:
    		"""
    		cd /data/MSc/2021/clipseq

    		STAR \\
        --runMode genomeGenerate \\
        --genomeDir index_star_ch20 \\
        --genomeFastaFiles /data/MSc/2021/clipseq/ref/chr20.fa \\
        --sjdbGTFfile /data/MSc/2021/clipseq/ref/chr20.gtf \\
        --genomeSAindexNbases 10 \\
        --outFileNamePrefix Hsapiens_chr20 \\
        --runThreadN 4
        """

}


////////////////////////////////////////////////////
/* --             CLIP PIPELINE                -- */
////////////////////////////////////////////////////

//Step 1 - FastQC

// Parameters

params.reads = "/data/MSc/2021/clipseq/sirna_trimmed_chr20.fq.gz"

// Initialise Channel
reads_ch = Channel.fromFilePairs(params.reads)

process FastQC {

      publishDir "${params.outdir}/QC/raw", mode:'copy'

      input:

      tuple val(key), file(reads) from reads_ch

      output:

      file("*.{html,zip}") into fastqc_ch

      script:
      """
      fastqc -q $reads
      """
}



//Step 2 - Trimming
// Test fasta file provided was aready trimmed

// Parameters

//params.reads = "/data/MSc/2021/clipseq/sirna_trimmed_chr20.fq.gz"

// Initialise Channel
//ch_reads= Channel.fromFilePairs(params.reads)

//process cutadapt {
    //tag "$key"
    //publishDir "${params.outdir}/cutadapt", mode: params.publish_dir_mode

    //input:
    //tuple val(key), file(reads) from ch_reads

    //output:
    //tuple val(key), path("${key}.trimmed.fq") into ch_trimmed
    //path "*.log" into ch_cutadapt_mqc

    //script:
    //"""
    //ln -s $reads ${key}.fq
    //cutadapt -j $task.cpus -a ${params.adapter} -m 12 -o ${key}.trimmed.fq ${key}.fq > ${key}_cutadapt.log
    //"""
//}


//Step 3 - Premapping


// Parameters



ch_trimmed = Channel.fromFilePairs(params.reads)


process premap {

    publishDir "${params.outdir}/premap", mode: 'copy'

    tag "$key"

    input:

    tuple val(key), file(reads) from ch_trimmed
    path(index) from ch_bt2_index.collect()

    output:
    tuple val(key), path("${key}.unmapped.fq.gz") into ch_unmapped
    tuple val(key), path("${key}.premapped.bam"), path("${key}.premapped.bam.bai")
    path "*.log" into ch_premap_mqc, ch_premap_qc

    script:
    """
    bowtie2 -p $task.cpus -x ${index[0].simpleName} --un-gz ${key}.unmapped.fq.gz -U $reads 2> ${key}.premap.log | \
    samtools sort -@ $task.cpus /dev/stdin > ${key}.premapped.bam && \
    samtools index -@ $task.cpus ${key}.premapped.bam
    """
}



//Step 4 - Aligning
process align {
    tag "$name"
    publishDir "${params.outdir}/mapped", mode: 'copy'

    input:
    tuple val(name), path(reads) from ch_unmapped
    path(index) from ch_star_index.collect()

    output:
    tuple val(name), path("${name}.Aligned.sortedByCoord.out.bam"), path("${name}.Aligned.sortedByCoord.out.bam.bai") into ch_aligned, ch_aligned_preseq
    path "*.Log.final.out" into ch_align_mqc, ch_align_qc

    script:
    clip_args = "--outFilterMultimapNmax 1 \
                --outFilterMultimapScoreRange 1 \
                --outSAMattributes All \
                --alignSJoverhangMin 8 \
                --alignSJDBoverhangMin 1 \
                --outFilterType BySJout \
                --alignIntronMin 20 \
                --alignIntronMax 1000000 \
                --outFilterScoreMin 10  \
                --alignEndsType Extend5pOfRead1 \
                --twopassMode Basic \
                --outSAMtype BAM Unsorted"
    """
    STAR \\
        --runThreadN $task.cpus \\
        --runMode alignReads \\
        --genomeDir $index \\
        --readFilesIn $reads --readFilesCommand gunzip -c \\
        --outFileNamePrefix ${name}. $clip_args

    samtools sort -@ $task.cpus -o ${name}.Aligned.sortedByCoord.out.bam ${name}.Aligned.out.bam
    samtools index -@ $task.cpus ${name}.Aligned.sortedByCoord.out.bam
    """
}


//Step 5 - Aligning QC

process preseq {
    tag "$name"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    input:
    tuple val(name), path(bam), path(bai) from ch_aligned_preseq

    output:
    path '*.ccurve.txt' into ch_preseq_mqc
    path '*.log'

    script:
    """
    preseq lc_extrap \\
        -output ${name}.ccurve.txt \\
        -verbose \\
        -bam \\
        -seed 42 \\
        $bam
    cp .command.err ${name}.command.log
    """
}


//Step 6 - Deduplicate AINE DID ALL OF THIS <3


process dedup {

        tag "$name"
        publishDir "${params.outdir}/dedup", mode: 'copy'

        input:
        tuple val(name), path(bam), path(bai) from ch_aligned

        output:
        tuple val(name), path("${name}.dedup.bam"), path("${name}.dedup.bam.bai") into ch_dedup, ch_dedup_pureclip, ch_dedup_rseqc
        //path "*.log" into ch_dedup_mqc, ch_dedup_qc

        script:
        """
        samtools rmdup -S $bam ${name}.dedup.bam
        samtools index -@ $task.cpus ${name}.dedup.bam
        """
}


//Step 7 - Identify crosslinks

process get_crosslinks {
    tag "$name"
    publishDir "${params.outdir}/xlinks", mode: 'copy'

    input:
    tuple val(name), path(bam), path(bai) from ch_dedup
    path(fai) from ch_fai_crosslinks.collect()

    output:
    tuple val(name), path("${name}.xl.bed.gz") into ch_xlinks_icount, ch_xlinks_paraclu, ch_xlinks_piranha
    tuple val(name), path("${name}.xl.bedgraph.gz") into ch_xlinks_bedgraphs
    path "*.xl.bed.gz" into ch_xlinks_qc

    script:
    """
    bedtools bamtobed -i $bam > dedup.bed
    bedtools shift -m 1 -p -1 -i dedup.bed -g $fai > shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' > pos.bed
    bedtools genomecov -dz -strand - -5 -i shifted.bed -g $fai | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' > neg.bed
    cat pos.bed neg.bed | sort -k1,1 -k2,2n | pigz > ${name}.xl.bed.gz
    zcat ${name}.xl.bed.gz | awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' | pigz > ${name}.xl.bedgraph.gz
    """
}

//Step 8d - Peak-call (Piranha)

process piranha_peak_call {

        tag "$name"
        publishDir "${params.outdir}/piranha", mode: 'copy'

        input:
        tuple val(name), path(xlinks) from ch_xlinks_piranha

        output:
        tuple val(name), path("${name}.${bin_size_both}nt_${cluster_dist}nt.peaks.bed.gz") into ch_peaks_piranha
        path "*.peaks.bed.gz" into ch_piranha_qc

        script:
        bin_size_both = 3
        cluster_dist = 3
        """
        pigz -d -c $xlinks | \\
        awk '{OFS="\t"}{for(i=0;i<\$5;i++) print }' \\
        > expanded.bed
        Piranha \\
            expanded.bed \\
            -s \\
            -b $bin_size_both \\
            -u $cluster_dist \\
            -o paraclu.bed
        awk '{OFS="\t"}{print \$1, \$2, \$3, ".", \$5, \$6}' paraclu.bed | \\
        pigz > ${name}.${bin_size_both}nt_${cluster_dist}nt.peaks.bed.gz
        """
}



//NOT SURE HOW TO DO DREME

//Step 8d2 - Motif (DREME)

//process piranha_motif_dreme {

        //tag "$name"
        //publishDir "${params.outdir}/piranha_motif", mode: 'copy'

        //input:
        //tuple val(name), path(peaks) from ch_peaks_piranha
        //path(fasta) from ch_fasta_dreme_piranha.collect()
        //path(fai) from ch_fai_piranha_motif.collect()

        //output:
        //tuple val(name), path("${name}_dreme/*") into ch_motif_dreme_piranha

        //script:
        //motif_sample = 1000
        //"""
        //pigz -d -c $peaks | awk '{OFS="\t"}{if(\$6 == "+") print \$1, \$2, \$2+1, \$4, \$5, \$6; else print \$1, \$3-1, \$3, \$4, \$5, \$6}' | \\
        //bedtools slop -s -l 20 -r 20 -i /dev/stdin -g $fai | \\
        //shuf -n $motif_sample > resized_peaks.bed
        //bedtools getfasta -fi $fasta -bed resized_peaks.bed -fo resized_peaks.fasta
        //dreme -norc -o ${name}_dreme -p resized_peaks.fasta
        //"""
//}

//Step 8 - QC plots

process clipqc {

    publishDir "${params.outdir}/clipqc", mode: 'copy'

    input:
    file ('premap/*') from ch_premap_qc.collect().ifEmpty([])
    file ('mapped/*') from ch_align_qc.collect().ifEmpty([])
    //file ('dedup/*') from ch_dedup_qc.collect().ifEmpty([])
    file ('xlinks/*') from ch_xlinks_qc.collect().ifEmpty([])
    file ('piranha/*') from ch_piranha_qc.collect().ifEmpty([])

    output:
    path "*.tsv" into ch_clipqc_mqc

    script:
    """
    /data/omccaffrey/MA5112/nextflow_assignment/running_clipseq.nf/clip_qc.py
    """
}


//Step 9 - MultiQC

process multiqc {

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    file (multiqc_config) from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc/*') from fastqc_ch.collect().ifEmpty([])
    file ('premap/*') from ch_premap_mqc.collect().ifEmpty([])
    file ('mapped/*') from ch_align_mqc.collect().ifEmpty([])
    path ('preseq/*') from ch_preseq_mqc.collect().ifEmpty([])
    file ('clipqc/*') from ch_clipqc_mqc.collect().ifEmpty([])

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"
    file "multiqc_plots"

    script:
    rtitle = ''
    rfilename = ''
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        rtitle = "--title \"${workflow.runName}\""
        rfilename = "--filename " + workflow.runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report"
    }
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $custom_config_file .
    """
}
