echo "chipmunk example for Sam (chip #3, Biotin Binders_1)"
# this determines where the output files go
OUTPUT=chipmunk_example_
APP=/home/dfeldman/s/app.sh
# if 1st and 2nd oligos are in separate files, interleave them before parsing
OLIGOS_A=/home/dfeldman/flycodes/assembly/chipmunk/chip-3-100000-230-pool-57.txt
OLIGOS_B=/home/dfeldman/flycodes/assembly/chipmunk/chip-3-100000-230-pool-58.txt
paste -d '\n' $OLIGOS_A $OLIGOS_B \
 | $APP parse-overlap-oligos stdin --name_col=None --dna_col=0 --output_prefix=$OUTPUT

echo "another example (Wei's CTLA4 library)"
# if 1st and 2nd oligos are in one file, use it directly
# name_col and dna_col correspond to columns in the oligos file
OLIGOS=/home/dfeldman/flycodes/assembly/CTLA4_L1.list
OUTPUT=Wei_example_
$APP parse-overlap-oligos $OLIGOS --name_col=0 --dna_col=1 --output_prefix=$OUTPUT
