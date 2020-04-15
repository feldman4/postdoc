reinitialize
initialize_settings

fetch 6D0TA
pdb BB1
align BB1, 6D0TA

hide (6D0TA)
hide (BB1)

show wire, BB1 and bb.
color gray, BB1

findpolar BB1 and bb., BB1_hbonds
color gray, BB1_hbonds

show wire, 6D0TA and bb.
chainbow 6D0TA
cnc 6D0TA

findpolar 6D0TA and bb., 6D0TA_hbonds

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

orient
deselect

_ set_view (\
_    -0.407253981,    0.899100959,   -0.160509929,\
_     0.566735625,    0.110962458,   -0.816391587,\
_    -0.716207802,   -0.423446298,   -0.554741442,\
_     0.000000000,    0.000000000, -118.366256714,\
_     6.546800613,    3.101176739,   38.186603546,\
_    51.626346588,  185.106201172,  -20.000000000 )

scene front, append

set_view (\
    -0.451695085,   -0.005935177,    0.892153859,\
    -0.004256834,    0.999978423,    0.004496935,\
    -0.892162204,   -0.001767335,   -0.451710612,\
     0.000000000,    0.000000000,  -81.323875427,\
     6.546800613,    3.101176739,   38.186603546,\
    14.583930969,  148.063827515,  -20.000000000 )

show sticks, 6D0TA_inner
show spheres, 6D0TA_inner
set sphere_transparency, 0.55, 6D0TA_inner

scene top, append

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

scene inside, append


#cmd.scene(key='front', action='recall', animate=0)
