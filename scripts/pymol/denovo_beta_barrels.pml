docstring = """
De novo beta barrel BB1 from Dou et al 2018. Designed from a schematic
representing strand pairing and glycine kinks. SSM data available.
"""

reinitialize
initialize_settings

fetch 6D0TA
pdb BB1.pdb
pdb BB1.tr.pdb
set_name BB1.tr, BB1_tr
align BB1, 6D0TA
align BB1_tr, 6D0TA

hide (BB1_tr)
hide (6D0TA)
hide (BB1)

show wire, BB1 and bb.
color gray, BB1

findpolar BB1 and bb., BB1_bb_hbonds
color gray, BB1_bb_hbonds

show wire, 6D0TA and bb.
chainbow 6D0TA
cnc 6D0TA

findpolar 6D0TA and bb., 6D0TA_bb_hbonds

select glycines, 6D0TA and (resi 26+44+54+56+80+104)
color yellow, glycines
show spheres, glycines and name CA
set sphere_scale, 0.3, glycines

findSurfaceAtoms 6D0TA
select exposed_atm_01, exposed_atm_01 and not bb.
select 6D0TA_inner, (6D0TA and polymer \
    and not (byres exposed_atm_01) and (not bb. or name CA))
delete exposed_atm_01

# matching residues, BB1 skips N-term methionine
stored.resi = []
iterate 6D0TA_inner, stored.resi.append(resi)
python
BB1_inner_resi = []
for i in set(stored.resi):
    BB1_inner_resi.append(str(int(i) - 1))
BB1_inner_resi = '+'.join(BB1_inner_resi)
python end

cmd.do(f'select BB1_inner, (BB1 and polymer and resi {BB1_inner_resi})')
cmd.do(f'select BB1_tr_inner, (BB1_tr and polymer and resi {BB1_inner_resi})')

# tryptophan corner
select 6D0TA_trp_corner, 6D0TA and resi 12+108
select BB1_trp_corner, BB1 and resi 11+107

orient
deselect

### FRONT VIEW

_ set_view (\
_    -0.407253981,    0.899100959,   -0.160509929,\
_     0.566735625,    0.110962458,   -0.816391587,\
_    -0.716207802,   -0.423446298,   -0.554741442,\
_     0.000000000,    0.000000000, -118.366256714,\
_     6.546800613,    3.101176739,   38.186603546,\
_    51.626346588,  185.106201172,  -20.000000000 )

scene front, append

### TOP VIEW

set_view (\
    -0.451695085,   -0.005935177,    0.892153859,\
    -0.004256834,    0.999978423,    0.004496935,\
    -0.892162204,   -0.001767335,   -0.451710612,\
     0.000000000,    0.000000000,  -81.323875427,\
     6.546800613,    3.101176739,   38.186603546,\
    14.583930969,  148.063827515,  -20.000000000 )

hide (BB1)
disable BB1_bb_hbonds

scene top, append

### CORNER

show sticks, 6D0TA_trp_corner
show wire, BB1_trp_corner

orient 6D0TA_trp_corner or BB1_trp_corner
label BB1_trp_corner and resi 107 and name CA, "R107"
label BB1_trp_corner and resi 11 and name CA, "W11"

scene corner, append

### CORE

label all, ""
show sticks, 6D0TA_inner
show spheres, 6D0TA_inner
set sphere_transparency, 0.15, 6D0TA_inner

orient 6D0TA_inner

scene core, append

### CORE PACKING

hide all
show wire, BB1_inner
show wire, 6D0TA_inner 
show mesh, BB1
set mesh_width, 0.1

nononpolarh

set_view (\
    -0.268817067,   -0.952537000,    0.142853782,\
     0.340915710,    0.044617876,    0.939034581,\
    -0.900840282,    0.301129431,    0.312744141,\
     0.000039103,    0.000007048, -110.444458008,\
     6.146798134,    2.708921432,   39.393108368,\
    92.897834778,  127.995033264,  -20.000000000 )

scene core_packing, append

### trRosetta

cmd.scene(key='front', action='recall', animate=0)

disable *hbonds
hide spheres
noh
as cartoon, BB1_tr
cartoon tube

color pink, BB1_tr

scene trRosetta, store

hide all, BB1_tr
as wire, BB1_tr and bb.

color gray, BB1
color green, 6D0TA
noh

scene trRosetta2, store

python

for i in range(1, 110):
    start = f'BB1_tr and resi {i} and name CA'
    end = f'6D0TA and resi {i + 1} and name CA'
    cmd.distance(f'dist_{i}', start, end)

cmd.group('distances', 'dist_*')
python end

hide labels
set dash_gap, 0.1, distances
set dash_length, 0.1, distances
set dash_width, 6, distances
disable BB1

#scene trRosetta3, store

hide everything, BB1_tr or 6D0TA or BB1
show cartoon, BB1_tr or 6D0TA or BB1
set cartoon_tube_radius, 0.1

disable BB1

label 6D0TA and name CA and resi 1, "N-term"
label 6D0TA and name CA and resi 111, "C-term"

scene trRosetta3, store
#cmd.scene(key='front', action='recall', animate=0)


python
def compare_from_to(start, end):
    cmd.do(f'create 6D0TA_cmp, 6D0TA and resi {start}-{end}')
    cmd.do(f'create BB1_tr_cmp, BB1_tr and resi {start - 1}-{end - 1}')
    compare_finish()

def compare_sel(first, second):
    cmd.do(f'create 6D0TA_cmp, {first}')
    cmd.do(f'create BB1_tr_cmp, {second}')
    compare_finish()

def compare_finish():
    cmd.do('hide all')
    cmd.do('align BB1_tr_cmp, 6D0TA_cmp')
    cmd.do('orient BB1_tr_cmp')
    cmd.do('as wire, 6D0TA_cmp')
    cmd.do('as wire, BB1_tr_cmp')
    cmd.do('deselect')
    cmd.do('cnc')
    cmd.do('noh')
    cmd.do('chainbow 6D0TA')
    cmd.do('group compare, BB1_tr_cmp 6D0TA_cmp 6D0TA_ref')

    cmd.do('chainbow 6D0TA')
    cmd.do('show cartoon, 6D0TA')

python end

#compare(30, 36)

compare_sel( \
    '( byres 6D0TA_inner  and x > 8)', \
    '( byres BB1_tr_inner and x > 8)', \
    )
#compare_sel('6D0TA_inner', 'BB1_tr_inner')
scene trRosetta4, store

print(docstring)