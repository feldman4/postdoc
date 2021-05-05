# sh akta_db_export_all_parallel.sh  2>&1 | tee export.log

JOBS=10

cd /net/scratch/dfeldman/akta_db/parallel
rm -rf all
mkdir all
cd all

AKTA_DB=/home/dfeldman/s/akta_db.sh

$AKTA_DB search Aaron
$AKTA_DB split_experiments chromatograms.csv
$AKTA_DB split_experiments /home/dfeldman/for/akta_db/chroma.hdf

touch commands.txt
work=`find export/ -name '*chromatograms.csv' | sort`
for f in $work; do
    echo $AKTA_DB export_and_plot $f >> commands.txt
done

parallel --jobs $JOBS < commands.txt

