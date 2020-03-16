log_open /Users/feldman/.pymol/log.pml,a

set cartoon_flat_sheets, 0
set cartoon_smooth_loops, 0

# run /path/to/home/pymol/load_sep.py

#   stick_radius -adjust thickness of atomic bonds
set stick_radius, 0.3

# save fetched PDB files here
set fetch_path, /tmp

fetch 6mrr
orient all


# util.chainbow("object-name")

python
def hide_water():
    cmd.hide('everything', 'resn hoh')

def chainbow(selection='all'):
    util.chainbow(selection)
python end

cmd.extend("nowater",hide_water)
cmd.extend("chainbow",chainbow)

nowater
chainbow