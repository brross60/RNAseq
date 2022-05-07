module load anaconda3
conda activate mapping
cd xxx


### TRIMMING ADAPTERS ###
fastq_filenames=($(ls *.fastq))

for i in "${fastq_filenames[@]}"
do
  i_basename=$(basename $i .fastq)
  echo "trimming adapters from $i_basename ..."
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 22 -o "$i_basename"_22bp.trim.fastq "$i_basename".fastq > "$i_basename"_22bp.cutadapt_log.txt
done

###############################################

### REAGULAR MAPPING WITH Bowtie2 ###
trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 8 -x /storage/home/hcoda1/9/bross60/p-mwhiteley3-0/rich_project_bio-whiteley/ref_genomes/MAB/MABATCC19977 -q -U "$i_basename".trim.fastq -S MMAS_NCBI_"$i_basename".sam 2>> MMAS_NCBI_"$i_basename".bowtie_output.txt
done

##########################################################

### REMOVING HUMAN (HUMANGRCh38)READS WITH Bowtie2 ###version 2.4.2
#bowtie2-build -f GRCh38_latest_genomic.fna HUMANGRCh38

trim_filenames=($(ls *trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i trim.fastq)
  echo "mapping to human $i_basename ..."
  bowtie2 --end-to-end -p 20 -x /storage/home/hcoda1/9/bross60/p-mwhiteley3-0/rich_project_bio-whiteley/ref_genomes/HUMAN/HUMANGRCh38 -q -U "$i_basename".fastq -S toHuman_"$i_basename".sam --un noHuman_"$i_basename".fastq 2>> Human_"$i_basename".bowtie_output.txt
done

### REMOVEING DECOY (decoyGLBR) READS WITH Bowtie2 ### version 2.4.2
#bowtie2-build -f RedonePangenome_Apr12_2021_Jan16_2022.fasta decoyGLBR

filenames=($(ls *.fastq))

for i in "${filenames[@]}"
do
  i_basename=$(basename $i .fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 20 -x decoyGLBR -q -U "$i_basename".fastq -S noHuman_toDecoyGLBR_"$i_basename".sam --un noHuman_noDecoyGLBR_"$i_basename".fastq 2>> decoyGLBR_"$i_basename".bowtie_output.txt
done

mkdir todecoy
mv noHuman_toDecoyGLBR*.sam todecoy

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
