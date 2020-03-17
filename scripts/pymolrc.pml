log_open /Users/feldman/.pymol/log.pml,a

run $PACKAGES/postdoc/scripts/commands.pml

window box, 100, 40, 1200, 800

set cartoon_flat_sheets, 0
set cartoon_smooth_loops, 0

set stick_radius, 0.3
set fetch_path, $TMPDIR


run $PACKAGES/postdoc/scripts/6mrr.pml