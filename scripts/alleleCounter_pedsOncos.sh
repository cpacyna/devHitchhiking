#!/bin/bash
###########
# script for standard alleleCounter stuff
##########

module load samtools/
module load alleleCount/

#bsub -env 'all, list='samples.txt'' -W 360 -G team154-vwork -o $PWD/AC.log.%J -e $PWD/AC.err.%J -R "span[hosts=1]" -R'select[mem>12000] rusage[mem=12000]' -M 12000 /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/genomic_analysis/unmatched/cgpvaf/alleleCounts/alleleCountSend.sh

for bampath in `ls /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/genomic_analysis/Bams/PD4*.bam`; do
  basename=`echo $(basename "$bampath") | cut -d'.' -f1`
  echo $basename
  bsub -env 'all' -G team154-vwork -J $basename -o $PWD/logs/${basename}.ac.log.%J -e $PWD/logs/${basename}.ac.err.%J -R "span[hosts=1]" -R'select[mem>12000] rusage[mem=12000]' -M 12000 /software/CASM/modules/installs/allelecount/shim-bin/alleleCounter -l /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/pediatric_mutations.txt -b $bampath -o /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/${basename%.*}_allelereport.txt -f 0 -F 0 -m 20
  #echo "/software/CASM/modules/installs/allelecount/shim-bin/alleleCounter -l /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/pediatric_mutations.txt -b $bampath -o /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/${basename%.*}_allelereport.txt -f 0 -F 0 -m 20"
done

# alleleCounter_pedsOncos.sh