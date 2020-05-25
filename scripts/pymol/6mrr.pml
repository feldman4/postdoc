reinitialize
initialize_settings

set sphere_scale, 0.2
set sphere_transparency, 0.4

fetch 6mrr

load $HOME/Downloads/supplementary_data/tested_foldit_design_pdbs/2002949_0000.pdb

select A1, /6mrr///14-27
select B1, /2002949_0000///11-24

align A1 and guide, B1 and guide

hide all
as wire, A1
as wire, B1
show spheres, (A1 or B1) and guide

nowater
noh

chainbow A1
color white, B1
cnc

orient A1

deselect