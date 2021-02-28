#!/bin/sh


# echo "Setting up"
# echo ./app.sh setup_files
# ./app.sh setup_files

# echo "Collecting input PDBs"
# echo ./app.sh load_monomer_table
# ./app.sh load_monomer_table


# LIMIT=300

# METHOD=sheets_only
# echo "Generating rpxdock commands for job submission"
# ./app.sh rpxdock_all --session $METHOD --method $METHOD --make_command \
#     --limit $LIMIT > rpxdock_commands_$METHOD.list

# METHOD=sheet_helix
# echo "Generating rpxdock commands for job submission"
# ./app.sh rpxdock_all --session $METHOD --method $METHOD --make_command \
#     --limit $LIMIT > rpxdock_commands_$METHOD.list


METHOD=sheets_only
METHOD=ss-EH_aa-ILVFM
LIMIT=30

for METHOD in ss-E_aa-ILVFM ss-EH_aa-ILVFM ss-EH_aa-ILFM ss-EH_aa-FM
do
    echo $METHOD
    ./app.sh rpxdock_all --session $METHOD --method $METHOD --limit $LIMIT \
        --make_command > sbatch/rpxdock_$METHOD.cmds

    # echo "Dumping symmetric pdbs"
    # ./app.sh dump_all_pdbs --session $METHOD --limit $LIMIT

    # echo "Exporting results"
    # ./app.sh export_results --session $METHOD

done