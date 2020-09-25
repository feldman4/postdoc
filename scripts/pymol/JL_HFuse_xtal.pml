# docstring = """ \
# Homo-trimer from Boyken et al 2016. Contains two copies of the same \
# hydrogen bond network. \
# """ 

reinitialize
initialize_settings

python
home = 'home/jlubner/working/scaffolds/C3'
id = '0046'

# relaxed xtal
name1 = f'HFuse_pH192_3_{id}'
f1 = f'{home}/{name1}/{name1}.pdb'
cmd.do(f'load {f1}')
cmd.do(f'rename {name1}, {id}_xtal_relaxed')

# xtal
name2 = f'HFuse_pH192_3_{id}_C2'
f2 = f'{home}/{name1}/model/xtal/{name1}_C2.pdb'
cmd.do(f'load {f2}')
cmd.do(f'rename {name1}_C2, {id}_xtal')
cmd.do(f'remove {id}_xtal and ((chain E) or (chain B) or (chain F))')

cmd.do('remove not polymer')

# design

f3 = f'{home}/{name1}/model/design_model/{name1}.pdb'
cmd.do(f'load {f3}')
cmd.do(f'rename {name1}, {id}_design')

constraint = 'resi 400-'
cmd.do(f'tmalign {constraint} and {id}_design, {constraint} and {id}_xtal')
cmd.do(f'tmalign {constraint} and {id}_xtal_relaxed, {constraint} and {id}_xtal')
cmd.do(f'disable {id}_xtal_relaxed')
cmd.do('orient')

python end

color tv_green, 0046_xtal
color violet, 0046_design
labeltermini

scene trimer, append

##############################

orient chain A and *design

label
hide cartoon, not chain A
chainbow
color gray90, *xtal

scene monomer_align_center, append

##############################

cmd.do(f'copy {id}_design_al_edge, {id}_design')
cmd.do(f'disable {id}_design')

constraint = 'chain A and resi 1-40'
cmd.do(f'tmalign {constraint} and {id}_design_al_edge, {constraint} and {id}_xtal')

orient chain A and *design

scene monomer_align_edge, append
