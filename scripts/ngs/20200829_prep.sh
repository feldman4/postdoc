HOME="NGS/20200829_DF_AB"
PAT="AGTCGC(.*)AAGA[CG][CG]"

cd $HOME

# mkdir fastq_T1
# cp fastq_11/T1*R1*gz fastq_T1
# gunzip fastq_T1/*gz

# mkdir match_T1
cd fastq_T1
for f in *fastq
    do
    rg -o $PAT $f --replace '$1' > ../match_T1/$f.match
done

# wc -l fastq_T1/*fastq > T1.wc