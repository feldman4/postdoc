reinitialize

# load
load for/YK/20201019_MiSeq_analysis/pymol/RA_05052_NTF2_mutant_counts.pdb
load for/YK/20201019_MiSeq_analysis/pymol/RA_05052_NTF2.pdb

set_name RA_05052_NTF2_mutant_counts, by_mutant_count
set_name RA_05052_NTF2, binder

copy binder_no_color, by_mutant_count
hide cartoon, by_mutant_count

# select
select binder, binder_no_color
select active_site, binder_no_color and resn KRR
select shell, byres binder_no_color within 3 of active_site
deselect

# representations
hide all
show cartoon, binder_no_color
justpolarh

show spheres, by_mutant_count and name CA
set sphere_scale, 0.3

# color
color bluewhite, target
color palecyan, binder_no_color
spectrum b, yellow_red, by_mutant_count and name CA, 0, 3

cnc

# misc
set label_size, 30
set label_color, salmon
show sticks, active_site

# scenes
disable by_*
enable by_mutant_count
scene by_ssm, append, , 0

disable by_mutant_count

python
order = 'PGFWYQNILVADREHKSTMC'
import glob
files = glob.glob('for/YK/20201019_MiSeq_analysis/pymol/RA_05052_NTF2_by_res_?.pdb')
files = sorted(files, key=lambda x: order.index(x[-5]))
for f in files:
    name = f.split('/')[-1].replace('.pdb', '')
    print('NAME', name)
    cmd.do(f'load {f}')
    cmd.do(f'show spheres, {name} and name CA')
    cmd.do(f'spectrum b, yellow_red, {name} and name CA')

    cmd.do(f'disable RA_05052_NTF2_by_res*')
    cmd.do(f'enable {name}')
    cmd.do(f'scene {name}, append, , 0')

python end
scene by_ssm, recall

## view
set_view (\
    -0.707935214,    0.350314379,    0.613276482,\
    -0.648484945,    0.021607704,   -0.760920405,\
    -0.279813111,   -0.936382592,    0.211877435,\
     0.000000000,    0.000000000, -111.755622864,\
    17.258131027,  -12.902084351,   -9.000738144,\
    82.358116150,  141.153137207,  -20.000000000 )
