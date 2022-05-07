module load anaconda3
conda activate mapping
cd xxx

#################################

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 20 -x decoyGLBR -q -U "$i_basename".trim.fastq -S decoyGLBR_"$i_basename".sam --un non-decoyGLBR_"$i_basename".fastq 2>> decoyGLBR_"$i_basename".bowtie_output.txt
done

###############################################

### REAGULAR MAPPING WITH Bowtie2 ###
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 8 -x /nv/hp10/bross60/data/ref_genomes/mmas.JM15300 -q -U "$i_basename".trim.fastq -S MMAS_NCBI_"$i_basename".sam 2>> MMAS_NCBI_"$i_basename".bowtie_output.txt
done

##########################################################

### REMOVING HUMAN (HUMANGRCh38)READS WITH Bowtie2 ###version 2.4.2
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end -p 20 -x /nv/hp10/bross60/data/ref_genomes/HUMAN/HUMANGRCh38 -q -U "$i_basename".trim.fastq -S toHuman_"$i_basename".sam --un noHuman_"$i_basename".fastq 2>> Human_"$i_basename".bowtie_output.txt
done

### REMOVEING DECOY (decoyGLBR) READS WITH Bowtie2 ### version 2.4.2
bowtie2-build -f RedonePangenome_Apr12_2021_Jan16_2022.fasta decoyGLBR

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 20 -x /nv/hp10/bross60/data/ref_genomes/DECOY/decoyGLBR -q -U noHuman_"$i_basename".fastq -S noHuman_toDecoyGLBR_"$i_basename".sam --un noHuman_noDecoyGLBR_"$i_basename".fastq 2>> decoyGLBR_"$i_basename".bowtie_output.txt
done

### MAPPING TO MAB (MABATCC19977) WITH Bowtie2 ### version 2.4.2
#bowtie2-build GCF_000195955.2_ASM19595v2_genomic.fna H37Rv

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 20 -x /nv/hp10/bross60/data/ref_genomes/MAB/MABATCC19977 -q -U noHuman_noDecoyGLBR_"$i_basename".fastq -S MAB_noHuman_noDecoy"$i_basename".sam 2>> MAB_noHuman_noDecoy_"$i_basename".bowtie_output.txt
done

##############################################################################

#OTHER REFS
#/nv/hp10/bross60/data/ref_genomes/TB/H37Rv

#################################################################################################################

### COUNTING FEATURES ###
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -g gene_id -t gene -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam
#-OR-
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -g gene_id -t CDS -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam


#change sam to bam so can view w/ ivg tools
samtools view -S -b input.sam > output.bam
samtools sort output.bam -o output_sort.bam
samtools index output_sort.bam
