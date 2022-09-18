mkdir dc_workshop
mkdir -p data/ref_genome
gunzip data/ref_genome/ecoli_rel606.fasta.gz
head data/ref_genome/ecoli_rel606.fasta
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
v sub/ ~/dc_workshop/data/trimmed_fastq_small
mkdir -p results/sam results/bam results/bcf results/vcf
bwa index data/ref_genome/ecoli_rel606.fasta
bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam
sudo apt-get install bcftools
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf --no-reference -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf
vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf
less -S results/vcf/SRR2584866_final_variants.vcf
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l
samtools index results/bam/SRR2584866.aligned.sorted.bam
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta


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
gunzip Ebere/dataset/references/hg19.chr5_12_17.fa.gz
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
