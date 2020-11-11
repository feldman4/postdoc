python
commands_pml = '$HOME/packages/postdoc/scripts/pymol/commands.pml'
cmd.extend('load_commands', lambda: cmd.run(commands_pml))
python end

load_commands

cmd.mouse('select_forward')
initialize_settings
