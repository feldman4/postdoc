# Find nearest matching designs for all assembled, translated protein sequences in NGS run.
# Runs in a couple minutes on one core.

QUERIES=/home/kipnis/match_and_design/xml_csts_cNTF2s/DNAworks/20201019_MiSeq_analysis/*.clean.txt
REF=/home/dfeldman/for/YK/20201019_MiSeq_analysis/designs.csv

# how to rename columns in the output table
# keep query info (i.e., read count)
QUERY_RENAME="1,count"
# keep reference info (i.e., design name)
REF_RENAME="dna,DROP,name,design_name"
# change match names from "reference" to "design"
RENAME="reference_match,design_match,reference_distance,design_distance,reference_equidistant,design_equidistant"

for QUERY in $QUERIES
do
    BASE=$(basename "$QUERY")
    echo "Processing $BASE"
    # col and header arguments used to parse input tables
    # keep arguments join results to input tables
    # rename argument renames or drops resulting columns
    APP=/home/dfeldman/s/app.sh
    $APP find_nearest_sequence $QUERY $REF \
        --col_query=0 \
        --col_reference=sequence --header_reference=True \
        --keep_query_table \
        --keep_reference_table \
        --rename_cols=$QUERY_RENAME,$REF_RENAME,$RENAME \
        > "$BASE.matched.csv"
done