delete all
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

select glycines, resn GLY and 6D0TA
color yellow, glycines
show spheres, glycines and name CA
set sphere_scale, 0.3

orient
deselect


_ set_view (\
_    -0.407253981,    0.899100959,   -0.160509929,\
_     0.566735625,    0.110962458,   -0.816391587,\
_    -0.716207802,   -0.423446298,   -0.554741442,\
_     0.000000000,    0.000000000, -118.366256714,\
_     6.546800613,    3.101176739,   38.186603546,\
_    51.626346588,  185.106201172,  -20.000000000 )
