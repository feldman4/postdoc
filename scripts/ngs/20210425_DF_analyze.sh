# rm -rf analysis
# mkdir analysis
# cp fastq/T1_B04_*R[12]_001.fastq.gz analysis
# cp fastq/T1_B05_*R[12]_001.fastq.gz analysis
# cp fastq/T1_B06_*R[12]_001.fastq.gz analysis
# gunzip analysis/*gz

rm -rf analysis_JL
mkdir analysis_JL

cp fastq/T1_A01_*R[12]_001.fastq.gz analysis_JL
cp fastq/T1_A02_*R[12]_001.fastq.gz analysis_JL
cp fastq/T1_A03_*R[12]_001.fastq.gz analysis_JL

gunzip analysis_JL/*gz
