fetch 7dfr, async=0
fetch 3dfr, async=0
nowater
align 3dfr, 7dfr
/cmd.set('grid_mode',1,'',0)
orient

cmd.select('sele',"byresi((((sele) or byresi((3dfr`1314))) and not ((byresi((3dfr`1314))) and byresi(sele))))",enable=1)
cmd.select('sele',"byresi((((sele) or byresi((3dfr`1364))) and not ((byresi((3dfr`1364))) and byresi(sele))))",enable=1)
set_name sele, ligands_3

deselect

cmd.select('sele',"byresi((((sele) or byresi((7dfr`1264))) and not ((byresi((7dfr`1264))) and byresi(sele))))",enable=1)
cmd.select('sele',"byresi((((sele) or byresi((7dfr`1237))) and not ((byresi((7dfr`1237))) and byresi(sele))))",enable=1)
set_name sele, ligands_7


deselect