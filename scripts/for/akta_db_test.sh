AKTA_DB=/home/dfeldman/s/akta_db.sh

rm -rf test
mkdir test
cd test

$AKTA_DB search koepnick
$AKTA_DB split_experiments chromatograms.csv

work=`find export/ -name '*chromatograms.csv' | sort`
for f in $work; do
    $AKTA_DB export_and_plot $f
done

