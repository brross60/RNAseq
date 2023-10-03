module load anaconda3
conda activate mapping



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
# -p 8 means use ad processors your PBS whould be at least nodes=1ppn=8 at minimum
# -very-sensistive is used for MAB but just means that the paramenter allow for better rmanning by seeding more. This sacrifies speed.
#-x is your reference genome indexes w/ the path if needed
# download your desired genome (.fna) then build your indexes w/ line below. You only need to do this once. This is what bowtie uses for mapping. Store them in your home directory not your scratch and use the path to call the data.
#2> outputs the bowtie stats in an txt file

#bowtie2-build -f MABATCC19977_genomic.fna MABATCC19977
#bowtie-inspect MABATCC19977 -s index-check.txt

trim_filenames=($(ls *.trim.fastq))

for i in "${trim_filenames[@]}"
do
  i_basename=$(basename $i .trim.fastq)
  echo "mapping $i_basename with bowtie2 ..."
  bowtie2 --end-to-end --very-sensitive -p 8 -x /storage/home/hcoda1/9/bross60/p-mwhiteley3-0/rich_project_bio-whiteley/ref_genomes/MAB/MABATCC19977 -q -U "$i_basename".trim.fastq -S MAB_NCBI_"$i_basename".sam 2>> MAB_NCBI_"$i_basename".bowtie_output.txt
done

### COUNTING FEATURES ###
# you much have featureCounts proggram in your folder and do >chmod +x featureCounts to make it executable
# -a is your reference (.gff) file downloaded from ncbi
#-t gene (inclides rRNA, tRNA, ncRNA)  or CDS only protein coding
#-g is what name to call the feature
#-O if the read maps to two genes/cds then thsi will count if for both. An operon will have RNA that map to muliple genes notre that only a single base overlap will make the read be counted to a gene.
# -s is the strandedness 0=unstranded 1=forward 2=reverse

./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -t gene -g gene_id -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam
#-OR-
./featureCounts -a /nv/hp10/bross60/data/ref_genomes/MAB/GCF_000069185.1_ASM6918v1_genomic.gff -t CDS -g gene_id -O -s 0 -o featureCounts_[insertname]_22bp.txt *.sam

### Calculates coverage with respect to MAB with samtools ###
#-@ use 16 threads for parallel processing
#-bS Convert SAM to BAM format
#-o output sorted BAM file
#-mA Skip alignments with MAPQ smaller than the specified minimum value (in this case, no minimum value is specified, so it'll default to 0)

for i in "${sam_filenames[@]}"
do
	i_basename=$(basename $i _22bp.sam)
	echo "analyzing $i_basename ..."
	samtools view -@ 16 -bS "$i_basename"_22bp.sam > "$i_basename"_22bp.bam
	samtools sort -@ 16 -o "$i_basename"_22bp.sorted.bam "$i_basename"_22bp.bam
	samtools coverage -mA "$i_basename"_22bp.sorted.bam
	samtools coverage -m -o "$i_basename"_22bp_coverage.txt "$i_basename"_22bp.sorted.bam
done






#################################################################################################################



#change sam to bam so can view w/ ivg tools
samtools view -S -b input.sam > output.bam
samtools sort output.bam -o output_sort.bam
samtools index output_sort.bam
