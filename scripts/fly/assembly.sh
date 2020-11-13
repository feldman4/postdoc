echo "CTLA4_L1"
app.sh parse-overlap-oligos CTLA4_L1.list
echo "FolateR1_L2"
app.sh parse-overlap-oligos FolateR1_L2.list --dna_col=0 --name_col=None
echo "MS library 1"
app.sh parse-overlap-oligos pool1.list --dna_col=0 --name_col=None
echo "MS library 2"
app.sh parse-overlap-oligos pool3.list