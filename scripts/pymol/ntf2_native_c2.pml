reinitialize

python

n = 1000
interface_cutoff = 3.9
exclude = '3er7',

# have to import consistently because namespace shared across scripts...
import glob

tmalign_exe = '/anaconda/bin/tmalign'

search = 'ntf2/symmetric_NTF2s_pdbs/*C2*pdb'
files = sorted(glob.glob(search))
files = [f for f in files if not any(x in f for x in exclude)]
names = [os.path.basename(f) for f in files]
print(f'Loading {len(files)} files from {search}')

def strip_common_prefix(xs):
    prefix = 0
    while all(x.startswith(xs[0][:prefix + 1]) for x in xs):
        prefix += 1
    return {x: x[prefix:] for x in xs}

def pymol_clean(name):
    if name.endswith('.pdb'):
        name = name[:-4]
    max_length = 14
    while name.startswith('_'):
        name = name[1:]
    return name[:max_length]

def show_interface(a, b, cutoff=interface_cutoff):
    name_a = a.replace(' ', '_')
    name_b = b.replace(' ', '_')
    int_a = name_a + '_int'
    int_b = name_b + '_int'
    limit = 'and not bb.'
    limit = ''
    cmd.do(f'create {int_a}, byres {a} {limit} within {cutoff} of {b} {limit}')
    cmd.do(f'create {int_b}, byres {b} {limit} within {cutoff} of {a} {limit}')
    cmd.do(f'show sticks, {int_a}')
    cmd.do(f'show sticks, {int_b}')
    #cmd.do(f'findpolar {int_a} to {int_b}')
    cmd.do(f'findpolar {int_a} to {int_b}, all')

    both = f'{int_a} or {int_b}'
    cmd.do(f'cnc {both}')
    cmd.do(f'set stick_radius, 0.1, ({both}) and bb.')
    cmd.do(f'set stick_radius, 0.25, ({both}) and (not bb. or name CA)')

cmd.do('set cartoon_gap_cutoff, 0')
cmd.do('set cartoon_rect_width, 0.1')
cmd.do('set cartoon_rect_length, 1')
cmd.do('set cartoon_oval_length, 1')
cmd.do('set ribbon_transparency, 0.7')
cmd.do('set cartoon_transparency, 0.5')
cmd.do('set sphere_transparency, 0.8')
cmd.do('set scene_animation_duration, 0.3')

fixed = strip_common_prefix(names)
loaded_names = []
for i in range(len(files[:n])):
    name = pymol_clean(fixed[names[i]]).replace('_C2', '')
    loaded_names += [name]
    cmd.load(files[i], name + '_')

    if i > 0:
         cmd.do(f'tmalign {name}_, {loaded_names[0]}_, exe={tmalign_exe}')



for name in loaded_names:

    cmd.do(f'select {name}_A, {name}_ and chain A')
    cmd.do(f'select {name}_B, {name}_ and chain B')
    
    cmd.do(f'as cartoon, {name}_')
    show_interface(name + '_A', name + '_B')
    
    cmd.do(f'color white, {name}_A')
    cmd.do(f'color slate, {name}_B')
    cmd.do(f'crt {name}*int')
    cmd.do('cnc')
    #cmd.do(f'create {name}_int_sph, {name}*int')
    #cmd.do(f'as spheres, {name}_int_sph')
    #cmd.do(f'crt {name}_int_sph')
    
    cmd.do(f'group {name}, {name}*')
    cmd.do(f'zoom *_, -6')
    cmd.do(f'scene {name}, store')
    cmd.do(f'disable {name}')

cmd.do(f'scene {loaded_names[0]}, recall')


cmd.do('deselect')


python end