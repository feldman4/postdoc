source activate /home/dfeldman/.conda/envs/df

HOME=/home/dfeldman/NGS/20201024_DF

mkdir -p $HOME/merged
mkdir -p $HOME/fastqc

pear --threads 4 \
 -f $HOME/fastq/T1_D12_S48_R1_001.fastq.gz \
 -r $HOME/fastq/T1_D12_S48_R2_001.fastq.gz \
 -o $HOME/merged/D12.fastq

fastqc $HOME/fastq/T1_D12*

mv $HOME/fastq/*fastqc*/ $HOME/fastqc
