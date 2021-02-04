set label_size, 22

set specular, 0.2
set cartoon_flat_sheets, 0
set cartoon_smooth_loops, 0

set stick_radius, 0.3

set scene_animation_duration, 0.75
set movie_fps, 15

set mesh_quality, 5
set min_mesh_spacing, 0.5

python

import socket
host = socket.gethostname()
if host == 'line':
	cmd.do('window box, 100, 40, 1720, 1000')
else:
	cmd.do('window box, 100, 40, 1200, 800')


cmd.do('set fetch_path, {}'.format(os.path.join(pdbs_dir, 'fetch')))

import os
log_file = os.path.join(os.environ['HOME'], '.pymol', 'log.pml')
cmd.do(f'log_open {log_file}, a')
python end

pml axes.py

