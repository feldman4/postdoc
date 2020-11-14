reinitialize
set scene_animation_duration, 0.4

rcsb 1mel
split_chains 1mel
delete 1mel
rename 1mel_L, target
rename 1mel_A, nanobody

deselect

color gray90, target
chainbow nanobody

showinterface target, nanobody

zoom inter*, 8

scene sticks, append

show spheres, interface_nanobody
show spheres, target
color gray90, target

chainbow nanobody

set_view (\
     0.636674941,   -0.721801400,   -0.271364421,\
    -0.552920818,   -0.182014793,   -0.813108981,\
     0.537513971,    0.667733788,   -0.514988124,\
     0.000000000,    0.000000000, -142.014251709,\
    21.279447556,   24.868373871,   13.027057648,\
   111.965118408,  172.063385010,  -20.000000000 )

scene spheres, append


fetch 1dpx
rename 1dpx, target_alone

align target_alone, target
as ribbon
color magenta, target_alone
set ribbon_transparency, 0.5, target_alone
set_view (\
     0.135097951,   -0.088357486,   -0.986879885,\
    -0.970479846,    0.189027011,   -0.149778202,\
     0.199778438,    0.977988362,   -0.060210720,\
    -0.000106089,    0.000025269,  -80.511558533,\
    13.389055252,   23.379587173,   13.502068520,\
    53.650672913,  107.375694275,  -20.000000000 )

scene induced_fit, append

scene sticks, recall, animate=0
