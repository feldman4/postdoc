python

import glob
import os
import inspect

home = os.path.join(os.environ['HOME'], 'packages', 'postdoc')
scripts_dir = os.path.join(home, 'scripts', 'pymol')
pdbs_dir = os.path.join(home, 'resources/pdbs')

exclude_pml = ['pymolrc.pml', 'commands.pml', 'settings.pml']


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

def find_polar(selection='all', mode='all', name=None):
    tmp = 'temp123'
    cmd.select(tmp, selection)

    if name is None:
        name = '{}_polar_conts'.format('_'.join(selection.split()))
    if mode == 'all':
        cmd.dist(
            name,
            f'({tmp}) and not (solvent)',
            f'({tmp}) and not (solvent)',
            quiet=1,mode=2,label=0,reset=1);
    if mode == 'nobb':
        cmd.dist(
            name,
            f'({tmp}) and not (solvent) and not bb.',
            f'({tmp}) and not (solvent)',
            quiet=1,mode=2,label=0,reset=1);

    cmd.enable(name)
    cmd.delete(tmp)

def hide_hydrogens(selection='all'):
    cmd.hide('({} and hydro)'.format(selection))

def color_not_carbon(selection='all'):
    util.cnc(selection);

def color_by_chain(selection='all'):
    util.color_chains(selection);

def run_script(name=None):
    files = list_scripts()
    if name is None:
        files = [os.path.relpath(f, scripts_dir) for f in files]
        pretty_print('Available scripts:', files)
        return
        
    script = fuzzy_match(files, name)
    if script:
        print(f'Running {script}')
        cmd.run(script)
        # print(f'Finished running {script}')

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

def load_pdb_grid(search, max_to_grid=20):
    try:
        from natsort import natsorted
        files = natsorted(glob.glob(search))
    except ModuleNotFoundError:
        files = sorted(glob.glob(search))

    cmd.do('set scene_animation_duration, 0')
    print(f'Loading {len(files)} files...')
    first_name = None
    for i, file in enumerate(files):
        name = os.path.basename(file)
        cmd.load(file, name)
        cmd.do(f'set grid_slot, {i+1}, {name}')
        if i % max_to_grid == max_to_grid - 1:
            cmd.do(f'scene {i}, store')
            cmd.do('disable all')
        if i == 0:
            first_name = name
        else:
            cmd.do(f'tmalign {name}, {first_name}')
    if i % max_to_grid != max_to_grid - 1:
        cmd.do(f'scene {i}, store')
        cmd.do('disable all')
    if i >= max_to_grid:
        cmd.do(f'scene {max_to_grid - 1}, recall')
    cmd.do('set grid_mode, 1')
    cmd.do('zoom, buffer=-10')

def color_by_residue_type(selection='all'):
    acid = 'oxygen'
    basic = 'neptunium'
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
        cmd.do(f'color {color}, backbone and resn {res}')

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

commands = [
# aliases
('rename', rename_selection),
# visibility
('nowater', hide_water),
('nohoh', hide_water),
('noh', hide_hydrogens),
('justpolarh', show_polar_h),
# coloring
('chainbow', chainbow),
('cnc', color_not_carbon),
('cbc', color_by_chain),
('crt', color_by_residue_type),
('glyballs', glycine_ca_spheres),
# selecting/labeling
('findpolar', find_polar),
('grabligands', select_ligands),
('labeltermini', label_termini),
# loading
('rcsb', fetch_with_defaults),
('pdbload', load_local_pdb),
('globload', load_pdb_grid),
# script management
('pmlrun', run_script),
('cml', list_commands),
('initialize_settings', initialize_settings),
('skeleton', skeleton),
]

for name, func in commands:
    cmd.extend(name, func)

load_external_scripts()

python end 
