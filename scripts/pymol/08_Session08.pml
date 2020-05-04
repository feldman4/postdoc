reinitialize
initialize_settings

python

working_dir = ('/Users/feldman/Downloads/Sessions/Session08_PackingAndRelax/'
    'export/')

cmd.load(working_dir + 'original.pdb', '2R0L')
cmd.load(working_dir + '2r0l_all_repack.pdb', 'all_repack')
cmd.load(working_dir + 'repack.pdb', 'H1_repack')
cmd.load(working_dir + 'redesign.pdb', 'redesign')


python end

#util.color_chains("(all and elem C)",_self=cmd)
orient

as ribbon
color gray, 2R0L
color red, all_repack
color green, H1_repack
color cyan, redesign
#cnc

python

template = 'select {0}_h1, {0} and chain H and res 27+28+29+30+37+38+39+40+41+42'
models = '2R0L', 'all_repack', 'H1_repack', 'redesign'
for model in models:
    cmd.do(template.format(model))

python end


#show wire, 2R0L_h1 and not elem h
#show wire, all_repack_h1 and not elem h
show sticks, redesign and not elem h
show sticks, H1_repack and not elem h

set stick_transparency, 0.5
set stick_radius, 0.2

set stick_transparency, 0, redesign
set stick_radius, 0.05, redesign

show spheres, H1_repack_h1 and name CA
color spheres, red
set sphere_scale, 0.3

# weird bug in pymol transparency, order matters?
disable H1_repack
enable H1_repack

deselect

orient H1_repack_h1
