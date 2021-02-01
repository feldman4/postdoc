reinitialize

python
import glob
import numpy as np

num_files = 30
max_name_length = 14
tmalign_exe = '/anaconda/bin/tmalign'

# names are too long
backbone_count = {}
def simplify_name(name):
    backbone = name.split('rd1')[0]
    backbone_count[backbone] = backbone_count.get(backbone, 0) + 1
    backbone_ix = list(backbone_count.keys()).index(backbone)
    sequence_ix = backbone_count[backbone]
    return f'{backbone_ix:03d}_{sequence_ix:03d}'

cmd.do('set scene_animation_duration, 0.3')
cmd.do('set cartoon_transparency, 0.5')

search = 'ntf2/201203_rocuronium_designs/*pdb'
exclude = []
files = sorted(glob.glob(search))
rs = np.random.RandomState(seed=0)
files = rs.choice(files, min(num_files, len(files)), replace=False)
files = sorted(files)

files = [f for f in files if not any(x in f for x in exclude)]
loaded_names = {f: simplify_name(f) for f in files}
first_name = list(loaded_names.values())[0]
print(f'Loading {len(files)} files from {search}')


for file, name in loaded_names.items():
    cmd.load(file, name + '_')
    if name != first_name:
         cmd.do(f'tmalign {name}_, {first_name}_, exe={tmalign_exe}')    

for name in loaded_names.values():
    cmd.do(f'v2 {name}')
    cmd.do(f'color bluewhite, {name} and name CA')
    cmd.do(f'group {name}, {name}*')
    cmd.do(f'zoom *_')
    cmd.do(f'scene {name}, store, view=0')
    cmd.do(f'disable {name}')

    cmd.do(f'set sphere_scale, 0.25, {name} and name CB')
    cmd.do(f'set stick_radius, 0.15, {name} and (name CA or name CB)')


cmd.do(f'scene {first_name}, recall')
cmd.do('deselect')

for file, name in loaded_names.items():
    print(' '*10, name, os.path.basename(file))

python end