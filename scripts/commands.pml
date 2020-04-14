python

import glob
import os
home = os.environ['HOME']
scripts_dir = 'drive/packages/postdoc/scripts'


def hide_water():
    cmd.hide('everything', 'resn hoh')

def chainbow(selection='all'):
    util.chainbow(selection)

def find_polar(selection='all'):
    tmp = 'temp123'
    cmd.select(tmp, selection)

    name = "{}_polar_conts".format(tmp)
    cmd.dist(
        name,
        "({}) and not (solvent)".format(tmp),
        "({}) and not (solvent)".format(tmp),
        quiet=1,mode=2,label=0,reset=1);
    cmd.enable(name)
    cmd.delete(tmp)

def hide_hydrogens(selection='all'):
    cmd.hide("({} and hydro)".format(selection))

def color_not_carbon(selection='all'):
    util.cnc(selection);

def run_script(name=None):
    if name is None:
        return list_scripts()
        
    if not os.path.exists(name):
        path = os.path.join(home, scripts_dir)
        files = glob.glob(os.path.join(path, '*pml'))
        files += glob.glob(os.path.join(path, '*py'))
        matches = [os.path.basename(f).startswith(name) for f in files]
        if sum(matches) > 1:
            print('Ambiguous name')
            return
        if sum(matches) == 0:
            print('No matching file')
            return
        name = files[matches.index(True)]

    cmd.run(name)

def list_scripts(absolute=False):
    path = os.path.join(home, scripts_dir, '*pml')
    files = glob.glob(path)
    exclude = ['pymolrc.pml', 'commands.pml']
    if not absolute:
        files = [os.path.basename(f) for f in files]
    print('*'*20)
    print('Available scripts:')
    for f in files:
        if any(f.endswith(x) for x in exclude):
            continue
        print('  ', f)
    print('*'*20)
    return files

def select_ligands(name='ligands'):
    selector = 'not pol. and not sol.'
    cmd.select(name, selector)

# local pdb loading
    
python end

cmd.extend("nowater", hide_water)
cmd.extend("nohoh", hide_water)
cmd.extend("noh", hide_hydrogens)
cmd.extend("chainbow", chainbow)
cmd.extend("findpolar", find_polar)
cmd.extend("cnc", color_not_carbon)
cmd.extend("pml_run", run_script)
cmd.extend("grabligands", select_ligands)