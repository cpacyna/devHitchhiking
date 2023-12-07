#!/bin/bash
###########
# script for standard alleleCounter stuff
##########

module load samtools/
module load alleleCount/

#bsub -env 'all, list='samples.txt'' -W 360 -G team154-vwork -o $PWD/AC.log.%J -e $PWD/AC.err.%J -R "span[hosts=1]" -R'select[mem>12000] rusage[mem=12000]' -M 12000 /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/genomic_analysis/unmatched/cgpvaf/alleleCounts/alleleCountSend.sh

for basename in `ls /lustre/scratch119/casm/team154pc/cp19/staged_adultOncos/Bams/2442`; do
  echo $basename
  bampath=/lustre/scratch119/casm/team154pc/cp19/staged_adultOncos/Bams/2442/$basename/mapped_sample/${basename}.sample.dupmarked.bam
  #bsub -env 'all' -G team154-vwork -J $basename -o $PWD/${basename}.ac.log.%J -e $PWD/${basename}.ac.err.%J -R "span[hosts=1]" -R'select[mem>12000] rusage[mem=12000]' -M 12000 /software/CASM/modules/installs/allelecount/shim-bin/alleleCounter -l /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/all_muts.txt -b $bampath -o /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/${basename%.*}_allelereport.txt -f 0 -F 0 -m 20
  command="/software/CASM/modules/installs/allelecount/shim-bin/alleleCounter -l /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/pediatric_mutations.txt -b $bampath -o /lustre/scratch119/casm/team274sb/cp19/projects/oncocytoma/DNA/alleleCounter/${basename%.*}_allelereport.txt -f 0 -F 0 -m 20"
  bsub -env 'all' -G team154-vwork -J $basename -o $PWD/${basename}.ac.log.%J -e $PWD/${basename}.ac.err.%J -R "span[hosts=1]" -R'select[mem>12000] rusage[mem=12000]' -M 12000 $command
done

# alleleCounter_adultOncos.sh