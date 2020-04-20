reinitialize
initialize_settings

fetch 5TJ3

nowater

show cartoon; show spheres, resn zn; color orange, resn tpo

select active, resn tpo
show sticks, active
orient active

select rosetta_hbonds, resi 100 or resi 162 or resi 164 or resi 486
show sticks, rosetta_hbonds

select polar_res, resn SER+THR+ASN+GLN+ARG+HIS+LYS+ASP+GLU


deselect
