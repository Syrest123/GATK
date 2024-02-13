#!/bin/bash
## ==Description==
This script downloads samples, performs quality control, indexes the ref, performs alignment; Provides the statistics for the alignment, sorts and indexes the bam, and performs variant calling and annotation
### Step 1: Dowloading samples

#SBATCH --job-name=GATK_SYRUS

#SBATCH --nodes=1
```bash
# Trial samples in samples.txt file
echo -e "SRR25434460\nSRR25434461" > triple.txt
sample=$(cat triple.txt)

for i in $sample
do
	echo "Downloading $i"
	fasterq-dump $i
done
```
### Step 2: Running FastQC & trimming
```bash
mkdir -p qualitycheck
qc_before="qualitycheck"
echo "Running FastQC..."

for i in $sample
do
	R1=${i}_R1.fastq.gz
	R2=${i}_R2.fastq.gz
#	fastqc $R1 $R2 -o $qc_before
	echo "fastqc for $i before trimming done"
#	fastp -i $R1 -I $R2 -o ${i}_trimmed_1.fastq -O ${i}_trimmed_2.fastq
	echo "Trimming for $i done"
#	fastqc ${i}_trimmed_1.fastq ${i}_trimmed_2.fastq
	echo "fastqc for $i after trimming done"
done
```
### Step 3(a): Downloading the reference genome
```bash
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta

#wget -nc https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.13.fa.gz
ref="Homo_sapiens_assembly38.fasta"
```
### 3(b): Indexing

```bash
# Indexing the reference sequence with bwa
echo "indexing reference"
mkdir -p reference
ref_dir="reference"
mv $ref $ref_dir
START="$( date +%s )"
if [ -e $ref_dir/*.ann ]
then
	echo " reference indexed"
else
	echo " Indexing reference"
	bwa index $ref_dir/$ref
	echo "Done indexing reference $ref"
fi
END="$( date +%s )"
echo "==========================Indexing the reference with bwa took: $[ $START - $END ] seconds==========================="

# Indexing reference with samtools to generate the .fai file

START="$( date +%s )"
if [ -e $ref_dir/*.fai ]
then
	echo "Reference indexed"
else
	echo "Indexing reference for gatk"
	samtools faidx $ref_dir/$ref
fi
END="$( date +%s )"
echo "==========================Indexing the reference with samtools took: $[ $START - $END ] seconds==========================="
```
### 3\(c): Creating .dicts file for the reference
```bash
START="$( date +%s )"
picard CreateSequenceDictionary --REFERENCE $ref_dir/$ref
END="$( date +%s )"
echo "==========================Generating dictionaries took: $[ $START - $END ] seconds=========================="
```
### Step 4(a): Downloading known sites
```bash
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget -nc https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz

# Unzipping them
gunzip 1000G_phase1.snps.high_confidence.hg38.vcf.gz
gunzip Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gunzip Homo_sapiens_assembly38.known_indels.vcf.gz

# Variable
known1="1000G_phase1.snps.high_confidence.hg38.vcf"
known2="Homo_sapiens_assembly38.known_indels.vcf"
known3="Mills_and_1000G_gold_standard.indels.hg38.vcf"
```
### 4(b): Indexing known sites
```bash
gatk IndexFeatureFile --input $known1
gatk IndexFeatureFile --input $known2
gatk IndexFeatureFile --input $known3
```

### Alignment, converting SAM to BAM, sorting BAM, index the sorted BAM and variant calling
```bash
# Setting up for alignment and calling variants
mkdir -p alignment_new
bam_dir="alignment_new"
mkdir -p gatk
gatk_dir="gatk"
```

### Step 5: Alignment step
```bash
for i in $sample
do
	if [ -e $bam_dir/${i}.bam ];
	then
		echo "----------------------Already aligned-----------------------"
	else
		START="$( date +%s )"
		R1=${i}_trimmed_1.fastq         # Assigning R1 to the corresponding read 1
		R2=${i}_trimmed_2.fastq         # Assigning R2 to the corresponding read 2
		bwa mem -aM -t 16 $ref_dir/$ref $R1 $R2 -R "@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA" | samtools view -Shb - > $bam_dir/${i}.bam          # Generating a bam file directly
		END="$( date +%s )"
		echo "==========================Alignment took: $[ $START - $END ] seconds==========================="
	fi
done
```

### Step 6: Sorting and indexing the bam
```bash
for i in $sample
do
	if [ -e $bam_dir/${i}_sorted.bam ] && [ -e $bam_dir/${i}_sorted.bam.bai ];
	then
		echo "-------------------Already sorted and indexed-------------------"
	else
		START="$( date +%s )"
	#samtools view -O BAM -o $dir1/${i}.bam $dir1/${i}.sam     # Converting the sam into bam
		samtools sort $bam_dir/${i}.bam > $bam_dir/${i}_sorted.bam      # Sorting the bam file
		samtools index $bam_dir/${i}_sorted.bam                     # Indexing the sorted bam
		END="$( date +%s )"
		echo "==========================Sorting and indexing $i took: $[ $START - $END ] seconds========================="
	fi
done
```

### Step 7: Generating aligment statistics
```bash
for i in $sample
do
	if [ -e $bam_dir/${i}.stats.txt ];
	then
	echo "-------------------Statistics already exits-----------------"
	else
		START="$( date +%s )"
		samtools flagstat $bam_dir/${i}.bam > $bam_dir/${i}.stats.txt   # Easy to understand stats for the alignment
		END="$( date +%s )"
		echo "==========================Generating alignment statistics took: $[ $START - $END ] seconds==========================="
	fi
done
```
## Running gatk 
### Step 8: Marking duplicates
```bash
for i in $sample
do
	if [ -e $gatk_dir/${i}_output.mkdup ];
	then
		echo "--------------------Marking duplicates already performed--------------------"
	else
		START="$( date +%s )"
		gatk MarkDuplicates --INPUT $bam_dir/${i}_sorted.bam --OUTPUT $gatk_dir/${i}_output.mkdup --METRICS_FILE $gatk_dir/${i}_output.metrics.txt --REMOVE_DUPLICATES false --CREATE_INDEX true # Generating the mkdup file and the metrics statistics
		END="$( date +%s )"
        	echo "==========================Marking duplicates took: $[ $START - $END ] seconds=================================="
	fi
done
```

### Step 9: Performing base recalibration
```bash
for i in $sample
do
	if [ -e $gatk_dir/${i}_recal_data.table ];
	then
		echo " --------------- Base recalibration already performed -------------------"
	else
		START="$( date +%s )"
#	gatk BuildBamIndex --INPUT $dir2/${i}_output.mkdup --OUTPUT $dir2/${i}_output.mkdup.bam.bai  # Indexing the mkdup file
		gatk BaseRecalibrator --input $gatk_dir/${i}_output.mkdup --output $gatk_dir/${i}_recal_data.table --reference $ref_dir/$ref --known-sites $known1 --known-sites $known2 --known-sites $known3     # Generating the recalibration scores
		END="$( date +%s )"
        	echo "==========================Generating alignment statistics took: $[ $START - $END ] seconds======================================"
	fi
done
```
### Step 10: Applying recalibration
```bash
for i in $sample
do
	START="$( date +%s )"
	gatk ApplyBQSR --bqsr-recal-file $gatk_dir/${i}_recal_data.table --input $bam_dir/${i}_sorted.bam --output $gatk_dir/${i}_BQSR.bam --reference $ref_dir/$ref --create-output-bam-index true # Applying the recalibration scores to the sorted bam
        END="$( date +%s )"
        echo "==========================Applying recalibration took: $[ $START - $END ] seconds======================================"
done
```
### Step 11: Calling variants
```bash
START="$( date +%s )"
#	gatk BuildBamIndex --INPUT $dir2/${i}_BQSR.bam --OUTPUT $dir2/${i}_BQSR.bam.bai    # Indexing the recalibrated bam
gatk Mutect2 --input $gatk_dir/SRR25434460_BQSR.bam --tumor-sample SRR25434460 --input $gatk_dir/SRR25434461_BQSR.bam --normal-sample SRR25434461 --output $gatk_dir/${i}.vcf --reference $ref_dir/$ref   # Variant calling
END="$( date +%s )"
echo "========================== Calling variants took: $[ $START - $END ] seconds======================================"
```

### Step 12: Annotation
```bash
START="$( date +%s )"
java -jar snpEff/snpEff.jar download hg38                                    # Downloading the annotation database
java -jar snpEff/snpEff.jar hg38 $gatk_dir/${i}.vcf > $gatk_dir/${i}_ann.vcf         # Annotating the vcf
END="$( date +%s )"
echo "==========================Annotation took: $[ $START - $END ] seconds======================================"
```

    rm $dir1/${i}.sam         # Removing the sam to free up space
# ==Done==