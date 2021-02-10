#!/bin/sh

SESSION=sheets_only
LIMIT=300

echo "Setting up"
echo ./app.sh setup_files
./app.sh setup_files

echo "Collecting input PDBs"
echo ./app.sh load_monomer_table
./app.sh load_monomer_table

echo "Generating rpxdock commands for job submission"
./app.sh rpxdock_all --session $SESSION --method sheets_only --make_command \
    --limit $LIMIT > rpxdock_commands.list

echo "Dumping symmetric pdbs for all results"
echo ./app.sh dump_all_pdbs --session $SESSION --limit_docks 5
./app.sh dump_all_pdbs --session $SESSION --limit_docks 5

