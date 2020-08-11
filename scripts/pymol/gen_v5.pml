python

n = 10

import glob
import os

experiments = {
'sample': 'red',
'sample_pssm': 'cyan',
'msa': 'yellow',
'vanilla': 'orange',
}
first_name = None
for exp in experiments:
	files = sorted(glob.glob(f'wfc/gen_v5/{exp}/*pdb'))
	files = [f for f in files if 'scwrl4' not in f]
	for file in files[:n]:
		name = os.path.basename(file)
		cmd.load(file, name)
		cmd.do(f'color {experiments[exp]}, {name}')
		if first_name is None:
			first_name = name
		else:
			cmd.do(f'align {name}, {first_name}')	

cmd.do('set grid_mode, 1')		
cmd.do('zoom, buffer=-10')

python end