#!/bin/sh

AKTA_DB=/home/dfeldman/s/akta_db.sh
AKTA_DB_LOCATION=/home/dfeldman/for/akta_db
EXPORT_LOCATION=/net/scratch/dfeldman/akta_db
# EXPORT_LOCATION=/home/dfeldman/misc/test

jobs=10
after="100 years ago"

print_usage() {
  printf "Usage: akta_db_export_all.sh --after=\"1 week ago\" --jobs <number of parallel jobs>"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n akta_db_export_all -o a:j: --long after:,jobs: -- "$@")
if [ "$?" != "0" ]; then
  print_usage
fi
eval set -- "$PARSED_ARGUMENTS"

while : 
do
  case "$1" in
    --after) after="$2" ; shift 2;;
    --jobs) jobs="$2" ; shift 2;;
    --) shift; break ;;
    *) print_usage ;;
  esac
done


cd $EXPORT_LOCATION
# rm -rf all
mkdir all
cd all

echo "Updating hdf tables"
$AKTA_DB export_hdf --path=$AKTA_DB_LOCATION

echo "Splitting experiments..."
$AKTA_DB split_experiments /home/dfeldman/for/akta_db/chroma.hdf --after="$after"

echo "Exporting UV data and overlay plots..."
work=`find export/ -name '*chromatograms.csv' | sort`
rm -rf commands.txt
for f in $work; do
    echo $AKTA_DB export_and_plot $f >> commands.txt
done

parallel --jobs $jobs < commands.txt
