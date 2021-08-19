reinitialize
set auto_zoom, off

set_view (\
     0.128595203,    0.861094296,   -0.491915703,\
     0.708403170,    0.267364174,    0.653207123,\
     0.693993509,   -0.432474077,   -0.575621068,\
     0.000000000,    0.000000000, -396.197235107,\
     0.432933807,  -28.855762482,   33.762413025,\
   326.631958008,  465.762512207,  -20.000000000 )

rcsb 2B5I
remove not pol.

split_chains
delete 2B5I
rename 2B5I_A, il2
rename 2B5I_B, beta
rename 2B5I_C, gamma
rename 2B5I_D, alpha

showinterface beta, gamma, name=beta_gamma
showinterface il2, beta, name=il2_beta
showinterface il2, gamma, name=il2_gamma

chainbow il2
color white, il2 and not bb.
cnc il2

as wire, il2_beta_beta
as wire, il2_gamma_gamma

set_view (\
     0.128595203,    0.861094296,   -0.491915703,\
     0.708403170,    0.267364174,    0.653207123,\
     0.693993509,   -0.432474077,   -0.575621068,\
     0.000000000,    0.000000000, -396.197235107,\
     0.432933807,  -28.855762482,   33.762413025,\
   326.631958008,  465.762512207,  -20.000000000 )

scene native, store, view=0

########################################

hide everything

rcsb 6DG5
remove not pol.
align 6DG5 and chain B, beta
split_chains 6DG5
delete 6DG5

rename 6DG5_A, neo2
rename 6DG5_B, beta2
rename 6DG5_C, gamma2

#showinterface beta, gamma, name=beta_gamma
showinterface neo2, beta2, name=neo2_beta
showinterface neo2, gamma2, name=neo2_gamma
as wire, neo2_beta_beta
as wire, neo2_gamma_gamma

chainbow neo2
color white, neo2 and not bb.
cnc neo2

set_view (\
     0.128595203,    0.861094296,   -0.491915703,\
     0.708403170,    0.267364174,    0.653207123,\
     0.693993509,   -0.432474077,   -0.575621068,\
     0.000000000,    0.000000000, -396.197235107,\
     0.432933807,  -28.855762482,   33.762413025,\
   326.631958008,  465.762512207,  -20.000000000 )

scene neo2, store, view=0
scene neo2_v, store, view=1

########################################

hide everything

rcsb 6DG6
rename 6DG6, neo2_alone
align neo2_alone, neo2

color gray, neo2_alone
set cartoon_transparency 0.5, neo2_alone

#v2 neo2
#v2 neo2_alone
show ribbon, neo2 or neo2_alone
hide cartoon
show wire, (name CA or not bb.) and (neo2 or neo2_alone)
cnc neo2 or neo2_alone

orient neo2
scene neo2_alone, store

# inaccurate loop
# loop KTTASEDEQEE
findseq KTTASEDEQEE, neo2_alone, loop_alone
findseq KTTASEDEQEE, neo2, loop_bound
deselect

align neo2 and loop_bound, neo2_alone and loop_alone
show spheres, (loop_alone or loop_bound) and name CA
set sphere_scale, 0.2
orient neo2
scene neo2_loop, store
