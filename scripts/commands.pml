python
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
    import os
    home = os.environ['HOME']
    scripts_dir = 'drive/packages/postdoc/scripts'
    path = os.path.join(home, scripts_dir, name)
    cmd.run(path)

python end

cmd.extend("nowater",hide_water)
cmd.extend("noh",hide_hydrogens)
cmd.extend("chainbow",chainbow)
cmd.extend("findpolar",find_polar)
cmd.extend("cnc", color_not_carbon)
cmd.extend("run_script", run_script)