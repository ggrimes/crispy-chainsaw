
/*
input params
*/

params.sra='SRR4204500'
params.kmer=21
params.outdir='results'
kmers = [19, 23, 25, 27, 29, 31]
//kmers = [19]

log.info """\
         R N A S E Q - N F   P I P E L I N E
         ===================================
         SRA:         : ${params.outdir}
         outdir:      : ${params.outdir}
         kmers        : ${kmers}
         """
         .stripIndent()

/*sra_ch = Channel
    .fromSRA(params.sra)
    .view()
*/

sra_ch = Channel
        .fromFilePairs("${baseDir}/data/raw_reads/SRR4204500/*{1,2}.fastq.gz")
        .view()



process cutadapt {

  publishDir "${params.outdir}/cutadapt"

  input:
  tuple(val(sampleID),path(fastqs)) from sra_ch

  output:
  tuple(val(sampleID),path("*cutadapt*fastq")) into trimmed_ch

  script:
  """
  cutadapt \
  -a AGATCGGAAGAGC \
  -A AGATCGGAAGAGC \
  -o ${sampleID}_cutadapt_1.fastq \
  -p ${sampleID}_cutadapt_2.fastq \
  ${fastqs}
  """
}

/*
velveth produces a hash table which is used by velvetg in the next step
velvetg builds and manipulates the De Bruijn graph and produces the assembly
*/
process velveth {

  publishDir "${params.outdir}/velveth"

  tag "velvet ${sampleID} ${kmer}"

  input:
  tuple(val(sampleID),path(fastq)) from trimmed_ch
  each(kmer) from kmers

  output:
  tuple(kmer,path("*contigs.fa")) into velveth_out_ch

  script:
  """
  velveth paired_${kmer} ${kmer} -shortPaired -fastq -separate  ${fastq}
  velvetg paired_${kmer}
  cp paired_${kmer}/contigs.fa ${sampleID}_${kmer}_contigs.fa
  """
}


process stats {

 publishDir "${params.outdir}/stats"

 input:
 tuple(kmer,path(contigs)) from velveth_out_ch

 output:
 path("${kmer}_stats.txt") into stats_ch

 script:
 """
 stats.sh ${contigs} > ${kmer}_stats.txt
 """
}


process combined_stats {

 publishDir "${params.outdir}/stats"

 input:
 path(stats) from stats_ch.collect()

 output:
 path("combined_stats.txt") into cout

 script:
 """
 grep -H 'contig N/L50' ${stats} > combined_stats.txt
 """
}
