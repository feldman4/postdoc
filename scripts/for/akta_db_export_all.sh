JOBS=10

cd /net/scratch/dfeldman/akta_db
rm -rf all
mkdir all
cd all

AKTA_DB=/home/dfeldman/s/akta_db.sh

echo "Splitting all experiments..."
$AKTA_DB split_experiments /home/dfeldman/for/akta_db/chroma.hdf

echo "Exporting UV data and overlay plots..."
work=`find export/ -name '*chromatograms.csv' | sort`
touch commands.txt
for f in $work; do
    echo $AKTA_DB export_and_plot $f >> commands.txt
done

parallel --jobs $JOBS < commands.txt
