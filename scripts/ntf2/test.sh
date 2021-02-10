#!/bin/sh

# use short name (sequence hash)
NAME_0=b507b9cdda
# use original name
NAME_1=hugh_ro9_5695_0_0.25_5_0_7404_0_0.35_1_W17_Y33_Q119d1m98a2_clean_0001_N22.rd1_000000232_0000400013_0000001_0.pdb
NAMES=$NAME,$NAME_1
SESSION=test

echo "Setting up"
echo ./app.sh setup_files
./app.sh setup_files

echo "Collecting input PDBs"
echo ./app.sh load_monomer_table
./app.sh load_monomer_table

echo "Docking one design locally"
echo ./app.sh rpxdock --session $SESSION --name $NAME_0 --method sheets_only
./app.sh rpxdock --session $SESSION --name $NAME_0 --method sheets_only

echo "Docking another design locally"
echo ./app.sh rpxdock --session $SESSION --name $NAME_1 --method sheets_only
./app.sh rpxdock --session $SESSION --name $NAME_1 --method sheets_only

echo "Generating rpxdock commands for job submission"
./app.sh rpxdock_all --session $SESSION --method sheets_only --make_command \
    --limit 100 > rpxdock_commands.list

echo "Dumping symmetric pdbs for one dock"
echo ./app.sh dump_pdbs --session $SESSION --name $NAME_0 --limit 5
./app.sh dump_pdbs --session $SESSION --name $NAME_0 --limit 5

echo "Dumping symmetric pdbs for all results"
echo ./app.sh dump_all_pdbs --session $SESSION --limit_docks 5
./app.sh dump_all_pdbs --session $SESSION --limit_docks 5

