case $1 in
    PDB30-20FEB17)
        # find /projects/ml/TrRosetta/PDB30-20FEB17/pdb -path '**/*.pdb'
        cat /projects/ml/TrRosetta/PDB30-20FEB17/list.small \
            | python -c "import sys; [print(f'/projects/ml/TrRosetta/PDB30-20FEB17/pdb/{line[1:3]}/{line[:-1]}.pdb') for line in sys.stdin]" \
            | pdb_stats.py \
            > rcsb/PDB30-20FEB17.small.csv
        ;;
    blast_cluster_30)
        ls rcsb/blast_cluster_30/*pdb | pdb_stats.py > rcsb/blast_cluster_30.csv
        ;;
    *)
        echo cases: PDB30-20FEB17, blast_cluster_30
        ;;
esac
