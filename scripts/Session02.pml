reinitialize
initialize_settings

fetch 5TJ3

nowater

show cartoon; show spheres, resn zn; color orange, resn tpo

select active, resn tpo
show sticks, active
orient * within 6 of active

select active_P, active and elem P
findpolar * within 4.5 of active_P, active_hbonds_pymol

select nearby_water, resn hoh within 4 of (active or elem Zn)
show spheres, nearby_water
set sphere_scale, 0.3, nearby_water
set sphere_transparency, 0.3, nearby_water
set sphere_color, cyan, nearby_water

cnc active

select rosetta_hbonds, resi 100 or resi 162 or resi 164 or resi 486
show sticks, rosetta_hbonds

_ set_view (\
_     0.617826164,   -0.694819748,   -0.368124932,\
_    -0.780047417,   -0.482594848,   -0.398281723,\
_     0.099079221,    0.533223391,   -0.840151906,\
_     0.000000000,    0.000000000,  -48.284255981,\
_    -3.611182928,  -29.466804504,   30.608322144,\
_    38.067676544,   58.500835419,  -20.000000000 )


scene active_site, append

select polar_residues, resn ARG+HIS+LYS+ASP+GLU+SER+THR+ASN+GLN
hide all
show spheres
color white
color orange, polar_residues
#set sphere_transparency, 0.5
set sphere_scale, 1

orient 5TJ3

scene balls, append

hide all
show sticks, bb.
orient 5TJ3

scene sticks, append

hide all
show cartoon

scene cartoon, append
deselect