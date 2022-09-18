mkdir -p Ebere/dataset
cd Ebere/dataset
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz
bash trim.sh
mv qc_reads trimmed-reads
repair.sh in1=trimmed_reads/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz in2=trimmed_reads/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz out1=SLGFSK-N_231335_r1_chr5_12_17_rep.fastq.gz out2=SLGFSK-N_231335_r2_chr5_12_17_rep.fastq.gz outsingle=single.fq
repair.sh in1=trimmed_reads/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz in2=trimmed_reads/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz out1=SLGFSK-T_231336_r1_chr5_12_17_rep.fastq.gz out2=SLGFSK-T_231336_r2_chr5_12_17_rep.fastq.gz outsingle=single.fq
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
mkdir -p results/sam results/bam results/bcf results/vcf
bwa index references/hg19.chr5_12_17.fa
bwa mem references/hg19.chr5_12_17.fa repaired/SLGFSK-N_231335_r1_chr5_12_17_rep.fastq.gz SLGFSK-N_231335_r2_chr5_12_17_rep.fastq.gz > results/sam/SLGFSK-N_231335.aligned.sam
bwa mem references/hg19.chr5_12_17.fa SLGFSK-T_231336_r1_chr5_12_17_rep.fastq.gz SLGFSK-T_231336_r2_chr5_12_17_rep.fastq.gz > results/sam/SLGFSK-T_231336.sam
samtools view -S -b results/sam/SLGFSK-N_231335.aligned.sam > results/bam/SLGFSK-N_231335.aligned.bam
samtools view -S -b results/sam/SLGFSK-T_231336.aligned.sam > results/bam/SLGFSK-T_231336.aligned.bam
samtools sort -o results/bam/SLGFSK-N_231335.aligned.sorted.bam results/bam/SLGFSK-N_231335.aligned.bam
samtools sort -o results/bam/SLGFSK-T_231336.aligned.sorted.bam results/bam/SLGFSK-T_231336.aligned.bam
samtools flagstat SLGFSK-N_231335.aligned.sorted.bam
samtools flagstat SLGFSK-T_231336.aligned.sorted.bam
bcftools mpileup -O b -o results/bcf/SLGFSK-N_231335_raw.bcf --no-reference -f references/hg19.chr5_12_17.fa results/bam/SLGFSK-N_231335.aligned.sorted.bam
bcftools mpileup -O b -o results/bcf/SLGFSK-T_231336_raw.bcf --no-reference -f references/hg19.chr5_12_17.fa results/bam/SLGFSK-T_231336.aligned.sorted.bam
bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-N_231335_variants.vcf results/bcf/SLGFSK-N_231335_raw.bcf
bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-T_231336_variants.vcf results/bcf/SLGFSK-T_231336_raw.bcf
vcfutils.pl varFilter results/vcf/SLGFSK-N_231335_variants.vcf  > results/vcf/SLGFSK-N_231335_final_variants.vcf
vcfutils.pl varFilter results/vcf/SLGFSK-T_231336_variants.vcf  > results/vcf/SLGFSK-T_231336_final_variants.vcf
less -S results/vcf/SLGFSK-N_231335_final_variants.vcf
less -S results/vcf/SLGFSK-T_231336_final_variants.vcf
grep -v "#" results/vcf/SLGFSK-N_231335_final_variants.vcf | wc -l
grep -v "#" results/vcf/SLGFSK-T_231336_final_variants.vcf | wc -l
