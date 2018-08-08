# /net/nfs/PAT/analysis/MPS-382/180601_BHV3L3BBXX/MedExome_hg19_Design_Annotation_files

# Input files
MedExome_hg19_empirical_targets.bed
MedExome_hg19_capture_targets.bed

bedtools="/net/nfs/PAT/home/matias/tools/bedtools2/bin/bedtools"

# FOR BOTH input files:

### 1 ### PRIMARY TARGETS

# Count the Number of Lines in the Input File
wc -l MedExome_hg19_empirical_targets.bed # 7720 194544

# Extract Region Coordinates (remove header line)
tail -n 194543 MedExome_hg19_empirical_targets.bed > MedExome_hg19_empirical_targets_coord_unmerged_unsorted.bed 

# Sort Region Coordinates
$bedtools sort -i MedExome_hg19_empirical_targets_coord_unmerged_unsorted.bed > MedExome_hg19_empirical_targets_coord_unmerged_sorted.bed

# Merge Overlapping and Book-Ended Region Coordinates
$bedtools merge -i MedExome_hg19_empirical_targets_coord_unmerged_sorted.bed > MedExome_hg19_primary_coord.bed


### 2 ### CAPTURE TARGETS
# Count the Number of Lines in the Input File
wc -l MedExome_hg19_capture_targets.bed # 223775

# Extract Region Coordinates (remove header line)
tail -n 223775 MedExome_hg19_capture_targets.bed > MedExome_hg19_capture_targets_unmerged_unsorted.bed 

# Sort Region Coordinates
$bedtools sort -i MedExome_hg19_capture_targets_unmerged_unsorted.bed > MedExome_hg19_capture_targets_unmerged_sorted.bed

#Merge Overlapping and Book-Ended Region Coordinates
$bedtools merge -i MedExome_hg19_capture_targets_unmerged_sorted.bed > MedExome_hg19_capture_coord.bed


################################################################
# 2. Create Picard Interval Lists
################################################################
# wd: /ccagc/home/matias/data/manifests/BCNHL_Seq_v2/makeManifest

###############
# 2A. Create Picard Target Interval Lists

module load samtools 

# take the header from a random bam file
SR50_head_sam='/net/nfs/PAT/analysis/MPS-382/CovMetrics_180601/OG_head.sam'

SR50_head_chr_bam='/net/nfs/PAT/analysis/MPS-382/CovMetrics_180601/LGG103_S1_chr.bam'
SR50_head_chr_sam='/net/nfs/PAT/analysis/MPS-382/CovMetrics_180601/SR50_head_chr_sam'
# and Create a Picard Interval List Header
samtools view -H /net/nfs/PAT/analysis/MPS-299/NHBCL_Seq_v1/160302_AC852VANXX/part1/bam/T15-18577_dedup_reads.bam > Random_bam_header.txt
samtools view -H $SR50_head_chr_bam > $SR50_head_chr_sam

#samtools view -H /net/nfs/PAT/analysis/MPS-299/NHBCL_Seq_v1/160302_AC852VANXX/part1/bam/T15-18577_dedup_reads.bam > Random_bam_header.txt

# create a picard interval list body
cat MedExome_hg19_primary_coord.bed | gawk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > MedExome_hg19_target_body.txt

# Concatenate to Create a Picard Target Interval List
#cat Random_bam_header.txt MedExome_hg19_target_body.txt > MedExome_hg19_target_intervals.txt

cat $SR50_head_sam MedExome_hg19_target_body.txt > MedExome_hg19_target_intervals.txt

sed -e 's/chr//g' MedExome_hg19_target_intervals.txt > MedExome_hg19_target_intervals_nochr.txt

cat $SR50_head_chr_sam MedExome_hg19_target_body.txt > MedExome_hg19_target_intervals_chr.txt


###############
# 2B. Create a Picard Bait Interval List

# Create a Picard Bait Interval List Body
cat MedExome_hg19_capture_coord.bed | gawk '{print $1 "\t" $2+1 "\t" $3 "\t+\tinterval_" NR}' > MedExome_hg19_bait_body.txt

# Concatenate to Create a Picard Bait Interval List
#cat Random_bam_header.txt MedExome_hg19_bait_body.txt > MedExome_hg19_bait_intervals.txt

cat $SR50_head_sam MedExome_hg19_bait_body.txt > MedExome_hg19_bait_intervals.txt

sed -e 's/chr//g' MedExome_hg19_bait_intervals.txt > MedExome_hg19_bait_intervals_nochr.txt

cat $SR50_head_chr_sam MedExome_hg19_bait_body.txt > MedExome_hg19_bait_intervals_chr.txt









