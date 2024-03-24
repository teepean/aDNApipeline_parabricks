#!/bin/bash

# Simple aDNA Pipeline for single-end fastqs

echo "*** Simple aDNA Pipeline for single-end fastqs ***"

check_program() {
  if ! command -v $1 &> /dev/null; then
    echo "$1 is required but not installed. Exiting."
    exit 1
  fi
}

# Check each required program
check_program samtools
check_program bwa
check_program fastp
check_program java

gpu_memory=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits --id=0 | awk '{print $1}')

# Check if the GPU memory is less than 12282 MiB
if [ "$gpu_memory" -lt 15000 ]; then
    echo "Not enough GPU memory"
    exit 1
fi

if [ "$#" -eq 0 ]; then
  echo "Syntax:"
  echo "    adnapipe_single <fastq1> <threads> <Population name> <Individual name>"
  echo
  exit 1
fi

if [ ! -e hs37d5.fa ]; then
  echo
  echo "You need to have hs37d5.fa in the main directory"
  echo "Do you want to download it? (y/n)"
  read -r choice
  case "$choice" in
    [Yy]* )
      echo "Downloading and uncompressing"
      curl --insecure -O https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
      gunzip -c hs37d5.fa.gz > hs37d5.fa
      rm hs37d5.fa.gz
      ;;
    [Nn]* ) 
      echo "Exiting"
      exit 1
      ;;
  esac
fi

if [ ! -e hs37d5.fa.fai ]; then
  echo
  echo "Indexing reference"
  samtools faidx hs37d5.fa
fi

if [ ! -e hs37d5.fa.bwt ]; then
  echo
  echo "You need to have hs37d5.fa bwa index in the main directory"
  echo "Do you want to index it? Notice: this might take several hours. (y/n)"
  read -r choice
  case "$choice" in
    [Yy]* )
      echo "Creating an index"
      bwa index hs37d5.fa
      ;;
    [Nn]* )
      echo "Exiting"
      exit 1
      ;;
  esac
fi

FASTQ=$1
THREADS=$2
POPNAME=$3
INDNAME=$4

if [ -d "$INDNAME" ]; then
  echo
  echo "Destination directory exists. Do you want to overwrite?"
  echo "This will delete the whole directory."
  read -p "Y/N " yn
  case $yn in
    [Yy]* )
      rm -rf "$INDNAME"
      mkdir -p "$INDNAME"
      ;;
    [Nn]* )
      echo
      echo "Quitting"
      echo
      exit
      ;;
    * ) echo "Please answer yes or no.";;
  esac
else
  mkdir -p "$INDNAME"
fi

echo "Removing adapters"

if [ -f "$PWD/$FASTQ" ]; then
  FASTQ1="$PWD/$FASTQ"
else
  FASTQ1=$FASTQ
fi

WORKPATH=$(dirname "$FASTQ1")

cd $INDNAME
fastp --in1 $FASTQ1 --out1  $INDNAME.fastp.fastq.gz --thread 4 --length_required 25 --json $INDNAME.fastp.json --html $INDNAME.fastp.html
cd ..

echo "Aligning"

#bwa aln -t $THREADS hs37d5.fa $INDNAME/$INDNAME.r1.fastq.gz -n 0.01 -l 1024 -k 2 > $INDNAME/$INDNAME.sai
#bwa samse -r "@RG\tID:ILLUMINA-$INDNAME\tSM:$INDNAME\tPL:illumina\tPU:ILLUMINA-$INDNAME-SE" hs37d5.fa $INDNAME/$INDNAME.sai $INDNAME/$INDNAME.r1.fastq.gz | samtools sort --no-PG -@ $THREADS -O bam - > $INDNAME/$INDNAME_SE.mapped.bam
docker run --gpus all --rm --volume "$WORKPATH":/workdir --volume $(pwd):/rootdir --volume $(pwd):/outputdir nvcr.io/nvidia/clara/clara-parabricks:4.2.1-1 pbrun fq2bam --low-memory --ref /rootdir/hs37d5.fa --in-se-fq /outputdir/${INDNAME}/${INDNAME}.fastp.fastq.gz --out-bam /outputdir/${INDNAME}/${INDNAME}_SE.mapped.bam
#samtools index -@ $THREADS ${INDNAME}/${INDNAME}_SE.mapped.bam

echo "Marking duplicates"

java -Xmx4g -jar picard/picard.jar MarkDuplicates INPUT=${INDNAME}/${INDNAME}_SE.mapped.bam OUTPUT=$INDNAME/$INDNAME_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE=$INDNAME/$INDNAME_rmdup.metrics VALIDATION_STRINGENCY=SILENT
samtools index -@ $THREADS $INDNAME/$INDNAME_rmdup.bam


echo "Damage profiling"
# Uncomment and adjust paths as necessary
# java -Xmx4g -jar bin/DamageProfiler-0.4.9.jar -i $INDNAME/$INDNAME_rmdup.bam -r hs37d5.fa -l 100 -t 15 -o . -yaxis_damageplot 0.30

echo "Trimming"
# Uncomment and adjust paths as necessary
# cygbin/bam trimBam $INDNAME/$INDNAME_rmdup.bam $INDNAME/tmp.bam -L 1 -R 1 
# bin/samtools sort -@ $THREADS $INDNAME/tmp.bam -o $INDNAME/$INDNAME.trimmed.bam 
# bin/samtools index -@ $THREADS $INDNAME/$INDNAME.trimmed.bam

echo "Genotyping"

samtools mpileup -B -q 30 -Q 30 -l v42.4.1240K.pos -f hs37d5.fa $INDNAME/$INDNAME_rmdup.bam | ./bin/pileupCaller --randomHaploid --sampleNames $INDNAME --samplePopName $POPNAME -f v42.4.1240K.snp -p $INDNAME/$INDNAME

echo "Done"
echo "The results are at $INDNAME subdirectory"
 
