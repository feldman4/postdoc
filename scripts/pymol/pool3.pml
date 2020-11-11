reinitialize

# load
load flycodes/pool3/pymol/by_ssm_entropy.pdb
load flycodes/pool3/pymol/by_res_P.pdb
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
show spheres, by_res_P and name CA
set sphere_scale, 0.3

# color
color bluewhite, target
color palecyan, binder_no_color
spectrum b, blue_white, by_ssm_entropy and name CA
spectrum b, yellow_red, by_pg_min and name CA
spectrum b, yellow_red, by_res_P and name CA

cnc

# misc
findpolar *interface, nobb
labeltermini binder
set label_size, 30
set label_color, salmon

# view
set_view (\
    -0.123192027,   -0.816734314,   -0.563709736,\
     0.937182486,   -0.282555521,    0.204572394,\
    -0.326360345,   -0.503098249,    0.800238013,\
     0.000000000,    0.000000000, -111.352996826,\
    -1.228801727,  -11.972537994,   -2.257324219,\
    87.791549683,  134.914443970,  -20.000000000 )
zoom *interface, 4


disable by_*
enable by_ssm_entropy
scene by_ssm, append, , 0

disable by_*
enable by_res_P
scene by_res_P, append, , 0
