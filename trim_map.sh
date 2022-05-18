module load anaconda3
conda activate mapping
#new pace updates keep througing me to the home folder. be sure to redirect just in case this happens.
cd xxx


### TRIMMING ADAPTERS ###
#This is a for loop to trim all your files which should look like *.fastq
#you do not need to change anything for this step
# -m 22 means keep only trimed reads that are >= 22bp 
# -o if name of output trim file "$i_basename" is a copy and paste of your orginal file name
# -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC this is the adaptor dont change

fastq_filenames=($(ls *.fastq))

for i in "${fastq_filenames[@]}"
do
  i_basename=$(basename $i .fastq)
  echo "trimming adapters from $i_basename ..."
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 22 -o "$i_basename"_22bp.trim.fastq "$i_basename".fastq > "$i_basename"_22bp.cutadapt_log.txt
done

###############################################
### REAGULAR MAPPING WITH Bowtie2 ###
# be sure to note your bowtie2 version!!!!!
# -p 8 means use ad processors your PBS whosul be nodes=1ppn=8 at minimum
# -very-sensistive is used for MAB but just means that teh paramenter allow for bette rmanning. this sacrifies speed
#-x is your reference genome indexes  w/ the path if needed
# download yor desired genome (.fna) then build your indexes w/ line below. You only need to do this once. this is what bowtie uses for mapping. 

#bowtie2-build -f MABATCC19977_genomic.fna MABATCC19977

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 8 -x /storage/home/hcoda1/9/bross60/p-mwhiteley3-0/rich_project_bio-whiteley/ref_genomes/MAB/MABATCC19977 -q -U "$i_basename".trim.fastq -S MMAS_NCBI_"$i_basename".sam 2>> MMAS_NCBI_"$i_basename".bowtie_output.txt
done

### COUNTING FEATURES ###
# you much gave featureCounts progrman in your filder and do chmod +x featureCounts to make it executable
# -a is your reference (.gff) file downloaded from ncbi
#-g is
#-t gene (inclides rRNA, tRNA, ncRNA)  or CDS only prootein coding
#-O of ythe read maps to two gens/cds then thsi will count if for both. An operon will have RNA that map to muliple genes
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -g gene_id -t gene -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam
#-OR-
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -g gene_id -t CDS -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam







#################################################################################################################



#change sam to bam so can view w/ ivg tools
samtools view -S -b input.sam > output.bam
samtools sort output.bam -o output_sort.bam
samtools index output_sort.bam
