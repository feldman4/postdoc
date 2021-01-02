HOME="NGS/20200829_DF_AB"
PAT_NTERM="AGTCGC(.*)AAGA[CG][CG]"
PAT_CTERM_DIALOUT="TCTCCCAAA(.*)"

cd $HOME

# mkdir fastq_T1
# cp fastq_11/T1*R1*gz fastq_T1
# gunzip fastq_T1/*gz

# mkdir match_T1
cd fastq_T1
for f in *fastq
    do
    rg -o $PAT_NTERM $f --replace '$1' > ../match_T1/$f.match_NTERM
done

for f in *fastq
    do
    rg -o $PAT_CTERM_DIALOUT $f --replace '$1' > ../match_T1/$f.match_CTERM
done


# wc -l fastq_T1/*fastq > T1.wc