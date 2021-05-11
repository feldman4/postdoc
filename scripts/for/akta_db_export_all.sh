#!/bin/sh

AKTA_DB=/home/dfeldman/s/akta_db.sh
AKTA_DB_LOCATION=/home/dfeldman/for/akta_db
EXPORT_LOCATION=/net/scratch/dfeldman/akta_db
# EXPORT_LOCATION=/home/dfeldman/misc/test
TIMESTAMP=`date +%s`

jobs=10
after="100 years ago"

print_usage() {
  printf "Usage: akta_db_export_all.sh --after=\"1 week ago\" --cores <number of parallel jobs>"
  exit 2
}

PARSED_ARGUMENTS=$(getopt -a -n akta_db_export_all -o a:j: --long after:,cores: -- "$@")
if [ "$?" != "0" ]; then
  print_usage
fi
eval set -- "$PARSED_ARGUMENTS"

while : 
do
  case "$1" in
    --after) after="$2" ; shift 2;;
    --cores) cores="$2" ; shift 2;;
    --) shift; break ;;
    *) print_usage ;;
  esac
done


cd $EXPORT_LOCATION
# rm -rf all
mkdir -p all
cd all

echo "Updating hdf tables"
$AKTA_DB export_hdf --path=$AKTA_DB_LOCATION

echo "Splitting experiments..."
$AKTA_DB split_experiments /home/dfeldman/for/akta_db/chroma.hdf --after="$after"



work=`find export/ -name '*chromatograms.csv' -mmin -10 | sort`
work_to_do=`echo "$work" | wc -l`
echo "Exporting UV data and overlay plots for $work_to_do experiments..."
commands="commands_$TIMESTAMP.txt"
for f in $work; do
    echo $AKTA_DB export_and_plot $f >> $commands
done

parallel --jobs $cores < $commands
