module load anaconda3
conda activate mapping




### MAPPING WITH Bowtie2 ###
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 8 -x /nv/hp10/bross60/data/ref_genomes/mmas.JM15300 -q -U "$i_basename".trim.fastq -S MMAS_NCBI_"$i_basename".sam 2>> MMAS_NCBI_"$i_basename".bowtie_output.txt
done

### COUNTING FEATURES ###
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MMAS_GCF_000497265.gff3 -g Parent -t CDS -s 1 -o featureCounts_MMAS_NCBI_EmoryHSP_18bp.txt MMAS_NCBI*.sam
Done



### MAPPING NON_DECOY TO MAB Bowtie2 ### version 2.4.2
bowtie2-build -f RedonePangenome_Apr12_2021_Jan16_2022.fasta decoyGLBR

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 20 -x decoyGLBR -q -U "$i_basename".trim.fastq -S decoyGLBR_"$i_basename".sam --un non-decoyGLBR_"$i_basename".fastq 2>> decoyGLBR_"$i_basename".bowtie_output.txt
done

### MAPPING WITH Bowtie2 ###
#bowtie2-build GCF_000195955.2_ASM19595v2_genomic.fna H37Rv
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 20 -x H37Rv -q -U non-decoyGLBR_"$i_basename".fastq -S TB_NCBI_"$i_basename".sam 2>> TB_H37Rv_"$i_basename".bowtie_output.txt
done





### MAPPING NON-HUMAN READS TO MAB WITH Bowtie2 ###
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 8 -x ref -q -U "$i_basename".trim.fastq -S Human_"$i_basename".sam --un Human_"$i_basename"_unmapped.fastq 2>> Human_"$i_basename".bowtie_output.txt
done


unmapped_filenames=($(ls *unmapped.fastq))

for i in "${unmapped_filenames[@]}"
do
  i_basename=$(basename $i .unmapped.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 8 -x ref -q -U Human_"$i_basename"22bp_unmapped.fastq -S MAB_NCBI_minus_MAC_"$i_basename".sam 2>> MAB_NCBI_minus_Human_"$i_basename".bowtie_output.txt
done









#change sam to bam so can view w/ ivg tools
samtools view -S -b input.sam > output.bam
samtools sort output.bam -o output_sort.bam
samtools index output_sort.bam
