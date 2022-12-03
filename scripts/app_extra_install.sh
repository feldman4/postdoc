# supplement /software/conda/envs/pyrosetta so app.sh can run
cd /home/dfeldman/packages/extra
pip install -t . natsort tqdm fire parse pyteomics lxml wget loess python-slugify dateparser
# these have lots of files and aren't needed
rm -r numpy* scipy* matplotlib* mpl*
