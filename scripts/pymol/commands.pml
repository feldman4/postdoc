python

import glob
import os
import inspect
from collections import defaultdict

home = os.path.join(os.environ['HOME'], 'packages', 'postdoc')
scripts_dir = os.path.join(home, 'scripts', 'pymol')
pdbs_dir = os.path.join(os.environ['HOME'], '.pymol/')
tmalign_exe = '/anaconda/bin/tmalign'

os.makedirs(pdbs_dir, exist_ok=True)

exclude_pml = ['pymolrc.pml', 'commands.pml', 'settings.pml']


def get_object_info(selection='all'):
    def accumulate(**kwargs):
        stored.arr.append(kwargs)

    stored.arr = []
    keys = 'model', 'chain', 'oneletter', 'resi', 'resv', 'resn', 'name', 
    kwargs = ', '.join([f'{k}={k}' for k in keys])
    space = {'accumulate': accumulate}
    cmd.iterate(f'({selection} & name CA)', f'accumulate({kwargs})', space=space)
    return stored.arr

def describe_chains(selection='all'):
    arr = get_object_info(selection)
    sequences = {}
    num_residues = {}
    
    grouped = defaultdict(list)
    for info in arr:
        grouped[(info['model'], info['chain'])] += [info]

    for model, chain in sorted(grouped):
        print(f'{model} -- chain {chain}')
        tmp = {x['resi']: x for x in grouped[(model, chain)]}
        by_res = sorted(tmp.values(), key=lambda x: x['resv'])
        min_resi, max_resi = min(tmp), max(tmp)
        selector = f'{model} & chain {chain} & polymer.protein'
        mass = util.compute_mass(selector, implicit=True)
        # not sure what the resi correspond to...
        #print(f'  {len(by_res)} residues (resi {min_resi} to {max_resi}), '
        #      f'MW {mass/1000:.2f} kDa')
        print(f'  {len(by_res)} residues, MW {mass/1000:.2f} kDa')
        print('  ' + ''.join(x['oneletter'] for x in by_res))
        
def list_pdbs():
    return glob.glob(os.path.join(pdbs_dir, '**/*pdb'), recursive=True)

def fuzzy_match(files, search):
    # could implement partial matching throughout directory structure
    # e.g., models/backbone1.pdb => mo/b
    files = sorted(set(files))

    # first exact matches
    if search in files:
        return search

    # now partial matches to filename
    matches = []
    for f in files:
        if os.path.basename(f).startswith(search):
            matches.append(f)

    if not matches:
        print('No matching file')
        return
    elif len(matches) > 1:
        print('Ambiguous:\n' + '\n'.join('  ' + x for x in matches))
        return
    else:
        return matches[0]

def load_local_pdb(name=None):
    files = list_pdbs()
    if name is None:
        files = [os.path.relpath(f, pdbs_dir) for f in files]
        pretty_print('Available pdbs:', files)
        return

    pdb = fuzzy_match(files, name)
    if pdb:
        print(f'Loading {pdb}')
        cmd.load(pdb)

def hide_water():
    cmd.hide('everything', 'resn hoh')

def chainbow(selection='all'):
    util.chainbow(selection)

def find_polar(selection='all', mode='nobb', name=None):

    split = selection.split(' to ')
    if len(split) == 2:
        selection_A, selection_B = split
    elif len(split) == 1:
        selection_A, selection_B = split[0], split[0]
    else:
        raise ValueError('more than two too many to split')

    tmp_A, tmp_B = 'temp123_A', 'temp123_B'
    cmd.select(tmp_A, selection_A)
    cmd.select(tmp_B, selection_B)

    if name is None:
        name = '{}_polar_conts'.format('_'.join(selection_A.split()))
    if mode == 'all':
        cmd.dist(
            name,
            f'({tmp_A}) and not (solvent)',
            f'({tmp_B}) and not (solvent)',
            quiet=1,mode=2,label=0,reset=1);
    elif mode == 'nobb':
        cmd.dist(
            name + 'A',
            f'({tmp_A}) and not (solvent) and not bb.',
            f'({tmp_B}) and not (solvent) and not bb.',
            quiet=1,mode=2,label=0,reset=1);
        cmd.dist(
            name + 'B',
            f'({tmp_A}) and not (solvent) and bb.',
            f'({tmp_B}) and not (solvent) and not bb.',
            quiet=1,mode=2,label=0,reset=1);
        cmd.dist(
            name + 'C',
            f'({tmp_A}) and not (solvent) and not bb.',
            f'({tmp_B}) and not (solvent) and bb.',
            quiet=1,mode=2,label=0,reset=1);
        # how to merge distance (?) objects?
        #cmd.do(f'create {name}, {name}*')
        #cmd.do(f'delete {name}A')
        #cmd.do(f'delete {name}B')
        #cmd.do(f'delete {name}C')
    elif mode == 'bbonly':
        cmd.dist(
            name,
            f'({tmp_A}) and not (solvent) and bb.',
            f'({tmp_B}) and not (solvent) and bb.',
            quiet=1,mode=2,label=0,reset=1);
        
    cmd.enable(name)
    cmd.delete(tmp_A)
    cmd.delete(tmp_B)

def hide_hydrogens(selection='all'):
    cmd.hide('({} and hydro)'.format(selection))

def color_not_carbon(selection='all'):
    util.cnc(selection);

def color_by_chain(selection='all'):
    util.color_chains(selection);

def color_yang(selection='all'):
    """Rosetta colors from Yang Hsia
    """
    colors = (
        ('teal', 'ASP+GLU'),
        ('teal', 'LYS+ARG'),
        ('teal', 'ASN+GLN'),
        ('teal', 'HIS'),
        ('marine', 'SER+THR'),
        ('grey90', 'ALA'),
        ('grey60', 'LEU+VAL+ILE+PHE'),
        ('grey60', 'MET'),
        ('yellow', 'CYS'),
        ('palegreen', 'TYR'),
        ('palegreen', 'TRP'),
        ('pink', 'GLY'),
        ('orange', 'PRO'),
        )

    for color, resn in colors:
        cmd.do(f'color {color}, {selection} and resn {resn}')


def run_script(name=None):
    if name == 'last' and hasattr(stored, 'last_script'):
        name = stored.last_script

    files = list_scripts()
    if name is None:
        files = [os.path.relpath(f, scripts_dir) for f in files]
        pretty_print('Available scripts:', files)
        return
        
    script = fuzzy_match(files, name)
    if script:
        print(f'Running {script}')
        stored.last_script = name
        cmd.run(script)


def list_scripts():
    files = glob.glob(os.path.join(scripts_dir, '*pml'))
    files += glob.glob(os.path.join(scripts_dir, '*py'))
    files += glob.glob(os.path.join(scripts_dir, 'external', '*py'))
    filtered = []
    for f in files:
        base = os.path.basename(f)
        if not any(base.startswith(x) for x in exclude_pml):
            filtered.append(f)
    return filtered

def pretty_print(header, items):
    print('*'*20)
    print(header)
    for x in items:
        print('  ', x)
    print('*'*20)

def select_ligands(name='ligands'):
    selector = 'not pol. and not sol.'
    cmd.select(name, selector)

def load_external_scripts():
    files = glob.glob(os.path.join(scripts_dir, 'external', '*py'))
    for f in files:
        cmd.do(f'run {f}')

def show_polar_h(selection='all', representation='wire'):
    """https://pymolwiki.org/index.php/Show
    """
    #cmd.do(f'hide everything, ele h and {selection}')
    #cmd.do(f'show lines, ele h and neighbor (ele n+o) and {selection}')
    # hide nonpolar hydrogens
    cmd.remove(f'{selection} & hydro')
    cmd.h_add(f'{selection} & (don.|acc.)')
    cmd.do(f'hide (h. and (e. c extend 1)) and {selection}')

def initialize_settings():
    cmd.run(os.path.join(scripts_dir, 'settings.pml'))

def label_termini(selection='all'):
    from collections import defaultdict
    stored.resi = defaultdict(list)
    cmd.iterate_state(-1, selection + ' and pol.', 'stored.resi[(model, chain)].append((int(resi), x))')
    for (model, chain) in stored.resi:
        first = min(stored.resi[(model, chain)])
        last = max(stored.resi[(model, chain)])
        select_this = f'{model} and chain {chain} and name CA'
        cmd.do(f'label {select_this} and resi {first[0]}, "N-term"')
        cmd.do(f'label {select_this} and resi {last[0]}, "C-term"')
        print(model, chain, first, last)

def list_commands():
    descriptions = []
    for name, f in commands:
        arguments = str(inspect.signature(f))
        descriptions += [name + arguments]
    pretty_print('Available commands:', descriptions)

def rename_selection(first, second=None):
    if second is None:
        old = 'sele'
        new = first
    else:
        old = first
        new = second

    cmd.do(f'set_name {old}, {new}')

def fetch_with_defaults(rcsb, assembly=1):
    cmd.do(f'fetch {rcsb}, type=pdb{assembly}')
    cmd.do(f'flatten_obj {rcsb}')
    cmd.do(f'delete {rcsb}')
    cmd.do(f'set_name flat, {rcsb}')
    hide_water()
    color_by_chain('not polymer.nucleic')

def load_pdb_grid(search, max_to_grid=20, max_to_load=None):
    try:
        from natsort import natsorted
        files = natsorted(glob.glob(search))
    except ModuleNotFoundError:
        files = sorted(glob.glob(search))
    if max_to_load:
        files = files[:int(max_to_load)]
    max_to_grid = int(max_to_grid)

    cmd.do('set scene_animation_duration, 0')
    print(f'Loading {len(files)} files...')
    names = [os.path.basename(f) for f in files]

    # remove shared prefixes
    prefix = 0
    if len(files) > 1:
        while len(set([x[:prefix + 1] for x in names])) == 1:
            prefix += 1

    first_name = None
    for i, file in enumerate(files):
        # random pymol limitations on object names
        # max length, can't start with underscore
        name = names[i][prefix:prefix+14]
        while name.startswith('_'):
            name = name[1:]

        cmd.load(file, name)
        cmd.do(f'set grid_slot, {i+1}, {name}')

        if i == 0:
            first_name = name
        else:
            cmd.do(f'tmalign {name}, {first_name}, exe={tmalign_exe}')
            cmd.do('orient')

        if i % max_to_grid == max_to_grid - 1:
            cmd.do(f'scene {i}, store')
            cmd.do('disable all')
    if i % max_to_grid != max_to_grid - 1 & i > 0:
        cmd.do(f'scene {i}, store')
        cmd.do('disable all')
    if i >= max_to_grid:
        cmd.do(f'scene {max_to_grid - 1}, recall')
    cmd.do('orient')
    cmd.do('zoom all, 4')
    cmd.do('set grid_mode, 1')
    rei;

def color_by_residue_type(selection='bb.'):
    acid = 'deepsalmon'
    basic = 'tantalum'
    nonpolar = 'vanadium'
    polar = 'thulium'
    cysteine = 'paleyellow'
    glycine = 'palecyan'

    colormap = {
        'asp': acid,
        'glu': acid,
        'arg': basic,
        'lys': basic,
        'his': basic,
        'met': nonpolar,
        'phe': nonpolar,
        'pro': nonpolar,
        'trp': nonpolar,
        'val': nonpolar,
        'leu': nonpolar,
        'ile': nonpolar,
        'ala': nonpolar,
        'ser': polar,
        'thr': polar,
        'asn': polar,
        'gln': polar,
        'tyr': polar,
        'cys': cysteine,
        'gly': glycine,
    }

    for res, color in colormap.items():
        cmd.do(f'color {color}, {selection} and resn {res}')

def glycine_ca_spheres(selection='all'):
    selector = f'name CA and {selection} and resn gly'
    cmd.do(f'show spheres, {selector}')
    cmd.do(f'set sphere_scale, 0.4, {selector}')
    # cmd.do(f'color palecyan, {selector}')
    
def skeleton(selection='all'):
    cmd.do(f'hide all, {selection}')
    cmd.do(f'show cartoon, {selection}')
    cmd.do(f'show ribbon, {selection}')
    cmd.do(f'show spheres, {selection} and name CA')
    cmd.do(f'show wire, {selection} and (name CA or name CB)')

    cmd.do(f'set line_width, 0.8, {selection}')
    cmd.do(f'set cartoon_transparency, 0.3, {selection}')
    cmd.do(f'set sphere_scale, 0.2, {selection}')

def axes_at_origin(scale=3):
    # axes.py
    from pymol.cgo import CYLINDER, cyl_text
    from pymol import cmd
    from pymol.vfont import plain

    # create the axes object, draw axes with cylinders coloured red, green,
    #blue for X, Y and Z
    s = scale
    obj = [
       CYLINDER, 0., 0., 0., s*10., 0., 0., s*0.2, 1.0, 1.0, 1.0, 1.0, 0.0, 0.,
       CYLINDER, 0., 0., 0., 0., s*10., 0., s*0.2, 1.0, 1.0, 1.0, 0., 1.0, 0.,
       CYLINDER, 0., 0., 0., 0., 0., s*10., s*0.2, 1.0, 1.0, 1.0, 0., 0.0, 1.0,
       ]

    # add labels to axes object (requires pymol version 0.8 or greater, I
    # believe

    #cyl_text(obj,plain,[-5.,-5.,-1],'Origin',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[s*10.,0.,0.],'X',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,s*10.,0.],'Y',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])
    cyl_text(obj,plain,[0.,0.,s*10.],'Z',0.20,axes=[[3,0,0],[0,3,0],[0,0,3]])

    # then we load it into PyMOL
    cmd.load_cgo(obj,'axes')

def show_interface(a, b, cutoff=3.5):
    name_a = a.replace(' ', '_')
    name_b = b.replace(' ', '_')
    cmd.do(f'select interface_{name_a}, byres {a} within {cutoff} of {b}')
    cmd.do(f'select interface_{name_b}, byres {b} within {cutoff} of {a}')
    cmd.do(f'show sticks, interface_{name_a}')
    cmd.do(f'show sticks, interface_{name_b}')
    cmd.do(f'findpolar interface_{name_a} to interface_{name_b}')
    both = f'interface_{name_a} or interface_{name_b}'
    cmd.do(f'cnc {both}')
    cmd.do(f'orient {both}')
    cmd.do(f'zoom {both}')
    cmd.do(f'set stick_radius, 0.15, ({both}) and bb.')
    cmd.do(f'set stick_radius, 0.25, ({both}) and (not bb. or name CA)')
    cmd.do(f'set cartoon_transparency, 0.8, {both}')


def show_stubs(selection='all', knobs=True):
    cmd.do(f'hide sticks, {selection}')
    cmd.do(f'show sticks, {selection} and name CA or name CB')
    if knobs:
        cmd.do(f'show spheres, {selection} and name CB')
        cmd.do(f'set sphere_scale, 0.3, {selection} and name CB')

def show_stubs_fancy(selection='all'):
    show_stubs(selection)
    color_by_residue_type(f'{selection} and name CB')

def external_wire(selection):
    gray = 'gray40'
    external = 'not (bb. or name CA or name CB)'
    cmd.do(f'show wire, {selection} and {external}')
    show_polar_h(f'{selection} and {external}')
    cmd.do(f'color {gray}, {selection} and {external}')
    color_not_carbon(f'{selection} and {external}')
    cmd.do(f'color {gray}, {selection} and {external} and elem H')


def view_one(selection='all'):
    cmd.do(f'hide everything, {selection}')
    cmd.do(f'show cartoon, {selection}')
    color_by_chain(selection)
    show_stubs_fancy(selection)


def view_two(selection='all'):
    cmd.do(f'hide everything, {selection}')
    cmd.do(f'show cartoon, {selection}')
    color_by_chain(selection)
    show_stubs_fancy(selection)
    external_wire(selection)
    find_polar(selection)


commands = [
# aliases
('rename', rename_selection),
# describe
('summarize', describe_chains),
# visibility
('nowater', hide_water),
('nohoh', hide_water),
('noh', hide_hydrogens),
('justpolarh', show_polar_h),
('show_stubs', show_stubs),
('stubs', show_stubs_fancy),
('v1', view_one),
('v2', view_two),
# coloring
('chainbow', chainbow),
('cnc', color_not_carbon),
('cbc', color_by_chain),
('crt', color_by_residue_type),
('cyang', color_yang),
('glyballs', glycine_ca_spheres),
# selecting/labeling
('findpolar', find_polar),
('grabligands', select_ligands),
('labeltermini', label_termini),
('showinterface', show_interface),
# loading
('rcsb', fetch_with_defaults),
('pdbload', load_local_pdb),
('globload', load_pdb_grid),
# script management
('pmlrun', run_script),
('cml', list_commands),
('initialize_settings', initialize_settings),
('skeleton', skeleton),
('labelorigin', axes_at_origin),
]

for name, func in commands:
    cmd.extend(name, func)

load_external_scripts()

python end 
