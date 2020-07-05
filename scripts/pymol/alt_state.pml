reinitialize
initialize_settings

load ~/from/CN/alt_state_pdbs/aligned_structures.pse

python
s = 'align {id}_min_e_low_rmsd_state///1-10/, {id}_alt_state///1-10/'

ids = [
 '997347_0007',
 '997523_0001',
 '2002689_0004',
 '2007855_0008',
 '2008226_0003',
]
for id in ids:
    cmd.do(s.format(id=id))

cmd.do('scene *, clear')
for id in ids:
    cmd.do('disable *')
    cmd.do(f'enable {id}_min_e_low_rmsd_state')
    cmd.do('orient')
    cmd.do('scene new, store')
    cmd.do('disable *')
    cmd.do(f'enable {id}_alt_state')
    cmd.do('scene new, store')

python end
