reinitialize
initialize_settings

docstring = """ \
Homo-trimer from Boyken et al 2016. Contains two copies of the same \
hydrogen bond network. \
""" 

rcsb 5J0H

show wire, not bb. or name CA or (resn PRO and name N)
show spheres, name CA
set sphere_scale, 0.18
set cartoon_transparency, 0.2
set line_width, 1

justpolarh
nohoh
cnc

findSurfaceAtoms
findpolar not exposed_atm_01, mode=nobb, name=5J0H_internal_polar
delete exposed_atm_01

python
positions = (
    'ASN`46',
    'LEU`60',
    'ASN`46',
    'ASN`45',
    'SER`29',
    'SER`22',
    'ASN`52',
    'ASN`53',
)
selector = ' or '.join(f'/5J0H///{x}' for x in positions)
cmd.do(f'show sticks, ({selector}) and (not bb. or name CA)')
python end
set stick_radius, 0.2


set_view (\
    0.998260856,    0.043887518,   -0.039328858,\
   -0.043827735,    0.999035418,    0.002387509,\
    0.039395750,   -0.000659464,    0.999222457,\
    0.000000000,    0.000000000,  -65.126068115,\
    0.000000000,   30.711643219,   98.851646423,\
   40.929851532,   89.322250366,  -20.000000000 )

scene 5J0H, append


deselect