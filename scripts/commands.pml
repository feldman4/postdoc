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
    name = "{}_polar_conts".format(selection)
    cmd.dist(
        name,
        "({}) and not (solvent)".format(selection),
        "({}) and not (solvent)".format(selection),
        quiet=1,mode=2,label=0,reset=1);
    cmd.enable(name)

def hide_hydrogens(selection='all'):
    cmd.hide("({} and hydro)".format(selection))

def color_not_carbon(selection='all'):
    util.cnc(selection);

def run_script(name):
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



python end

cmd.extend("nowater",hide_water)
cmd.extend("noh",hide_hydrogens)
cmd.extend("chainbow",chainbow)
cmd.extend("findpolar",find_polar)
cmd.extend("cnc", color_not_carbon)
cmd.extend("pml_run", run_script)
cmd.extend("pml_list", list_scripts)