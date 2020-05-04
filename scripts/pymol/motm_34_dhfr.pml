delete all
fetch 7dfr, async=0
nowater
grabligands
select interface, (byres pol. within 4 of ligands) and (not (backbone and not name CA))
select other_crap, pol. and not interface
color orange, other_crap
color cyan, interface

color green, ligands
cnc ligands

cmd.select('sele','none')
cmd.select('sele',"byresi((((sele) or byresi((7dfr`1256))) and not ((byresi((7dfr`1256))) and byresi(sele))))",enable=1)
util.cba(154,"sele",_self=cmd)

show sticks, interface
cnc interface

show surface, pol.
hide cartoon, pol.

findpolar interface or ligands

deselect

orient interface

set transparency, 0.3

_ set_view (\
_    -0.660127938,    0.118051484,    0.741818488,\
_     0.750476062,    0.145590603,    0.644661248,\
_    -0.031898290,    0.982275724,   -0.184705645,\
_     0.000000000,    0.000000000, -157.866043091,\
_    18.009531021,   23.188898087,   32.594112396,\
_  -513.267272949,  828.999328613,  -20.000000000 )
