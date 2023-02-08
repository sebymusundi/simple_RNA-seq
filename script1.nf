#! /usr/bin/env nextflow 

 
params.reads= "/home/seby/Afribop/Genomics/nf-phdcourse22/data/ggal/*_{1,2}.fq"
params.transcriptome_file ="/home/seby/Afribop/Genomics/nf-phdcourse22/data/ggal/transcriptome.fa"
params.annotation_file="/home/seby/Afribop/Genomics/nf-phdcourse22/data/ggal/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.108.chr.gtf"
params.outdir ="/home/seby/nextflow.files/results"


println "reads: $params.reads"
println "transcriptome: $params.transcriptome_file"
println "outdir: $params.outdir"
println "gtf_file: $params.annotation_file"

log.info """\
		RNA-SEQ PIPELINE
		---------------------------------------------
		transcriptome:  "$params.transcriptome_file"
		reads:	        "$params.reads"
		gtf_file:	"$params.annotation_file"
		outdir:		"$params.outdir"
		---------------------------------------------
		"""
		.stripIndent()
	

	
	
process FASTQC {
		conda '/home/seby/anaconda3/envs/bioinfo'
		publishDir "${params.outdir}"	
		
		input:
		tuple val(x), path(reads)

		output:
		path 'fastqc_out'

		script:
		"""
		mkdir fastqc_out
		fastqc -o fastqc_out $reads
		"""
}



process INDEX_HISAT{
		conda '/home/seby/anaconda3/envs/bioinfo'
		 publishDir "${params.outdir}"		
	
		input:
		path transcriptome

		output:
		tuple val("hisat_output"), path("hisat_output*")

		script:

		"""
		hisat2-build $transcriptome  hisat_output
		"""
}		

process ALIGN_HISAT{
		conda '/home/seby/anaconda3/envs/bioinfo'
		 publishDir "${params.outdir}"		

		input:
		tuple val(x), path(reads) , val(index), path(align_input)		
		
		output:
		tuple val(x), path("${x}.sam")

		script:
		"""
		hisat2 -fx $index -q  -1 ${reads[0]}   -2 ${reads[1]}   -S ${x}.sam
		"""
}
process SAMTOOLS_VIEW{
		conda '/home/seby/anaconda3/envs/bioinfo'
		 publishDir "${params.outdir}"		

		input:
		tuple val(x), path(sam_files)
				
		output:
		tuple val(x), path("${x}.sorted.bam.bai") , path("${x}.sorted.bam")

		script:
		"""
		samtools view -O BAM ${sam_files} -o ${x}.bam
		samtools sort ${x}.bam -o ${x}.sorted.bam
		samtools index ${x}.sorted.bam  
		"""
}
process READS_COUNT{
		 publishDir "${params.outdir}"

		input:
		tuple  path(gtf_file), val(x), path(bam_files_bai), path(bam_files)
		
		output: 
		path "*.tsv"

		script:
		"""
		Rscript /home/seby/Afribop/Genomics/nf-phdcourse22/data/ggal/subread_script.R ${bam_files} ${gtf_file} 
		"""
}


workflow{
my_reads=Channel.fromFilePairs("$params.reads")
//my_reads.view()
my_transcriptome=Channel.fromPath("$params.transcriptome_file")
my_gtf=Channel.fromPath("$params.annotation_file")
fastqc_ch=FASTQC(my_reads)
index_hisat_ch=INDEX_HISAT(my_transcriptome)
//index_hisat_ch.view()
//my_reads.combine(index_hisat_ch).view()
align_hisat_ch=ALIGN_HISAT(my_reads.combine(index_hisat_ch))
samtools_view_ch=SAMTOOLS_VIEW(align_hisat_ch)
//my_gtf.combine(samtools_view_ch).view()
READS_COUNT(my_gtf.combine(samtools_view_ch))
}
