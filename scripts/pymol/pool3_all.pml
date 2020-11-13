reinitialize

# load
load flycodes/pool3/pymol/by_ssm_entropy.pdb
load flycodes/pool3/pymol/target.pdb
copy binder_no_color, by_ssm_entropy
hide cartoon, by_ssm_entropy

# select
select binder, binder_no_color
select binder_interface, byres binder within 3 of target, 
select target_interface, byres target within 3 of binder
deselect

# representations
hide all
show cartoon, target or binder_no_color
show sticks, *interface and (not bb. or name CA)
justpolarh

show spheres, by_ssm_entropy and name CA
set sphere_scale, 0.3

# color
color bluewhite, target
color palecyan, binder_no_color
spectrum b, blue_white, by_ssm_entropy and name CA

cnc

# misc
findpolar *interface, nobb
labeltermini binder
set label_size, 30
set label_color, salmon

# scenes
disable by_*
enable by_ssm_entropy
scene by_ssm, append, , 0

python
order = 'PGFWYQNILVADREHKSTMC'
import glob
files = glob.glob('flycodes/pool3/pymol/by_res_?.pdb')
files = sorted(files, key=lambda x: order.index(x[-5]))
for f in files:
    name = f.split('/')[-1].replace('.pdb', '')
    print(name)
    cmd.do(f'load {f}')
    cmd.do(f'show spheres, {name} and name CA')
    cmd.do(f'spectrum b, yellow_red, {name} and name CA')

    cmd.do(f'disable by_*')
    cmd.do(f'enable {name}')
    cmd.do(f'scene {name}, append, , 0')

python end
scene by_ssm, recall

# view
set_view (\
    -0.123192027,   -0.816734314,   -0.563709736,\
     0.937182486,   -0.282555521,    0.204572394,\
    -0.326360345,   -0.503098249,    0.800238013,\
     0.000020284,    0.000000528, -118.481704712,\
    -2.449506283,  -13.798474312,   -5.327538490,\
    93.411872864,  143.551513672,  -20.000000000 )
zoom *interface, 4