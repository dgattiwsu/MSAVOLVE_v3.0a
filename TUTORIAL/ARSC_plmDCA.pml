
cmd.delete("all")
load ARSC_plmDCA.pdb
cmd.spectrum("b",selection=("ARSC_plmDCA"),quiet=0)
dist ////74/CA/,////82/CA/
select sele, resi 74+82
cmd.show("sticks","sele")
dist ////23/CA/,////25/CA/
select sele, resi 23+25
cmd.show("sticks","sele")
dist ////61/CA/,////62/CA/
select sele, resi 61+62
cmd.show("sticks","sele")
dist ////39/CA/,////81/CA/
select sele, resi 39+81
cmd.show("sticks","sele")
dist ////49/CA/,////94/CA/
select sele, resi 49+94
cmd.show("sticks","sele")
dist ////83/CA/,////90/CA/
select sele, resi 83+90
cmd.show("sticks","sele")
dist ////59/CA/,////65/CA/
select sele, resi 59+65
cmd.show("sticks","sele")
dist ////48/CA/,////101/CA/
select sele, resi 48+101
cmd.show("sticks","sele")
dist ////49/CA/,////51/CA/
select sele, resi 49+51
cmd.show("sticks","sele")
dist ////18/CA/,////19/CA/
select sele, resi 18+19
cmd.show("sticks","sele")
dist ////51/CA/,////56/CA/
select sele, resi 51+56
cmd.show("sticks","sele")
dist ////44/CA/,////48/CA/
select sele, resi 44+48
cmd.show("sticks","sele")
dist ////21/CA/,////26/CA/
select sele, resi 21+26
cmd.show("sticks","sele")
dist ////31/CA/,////34/CA/
select sele, resi 31+34
cmd.show("sticks","sele")
dist ////53/CA/,////76/CA/
select sele, resi 53+76
cmd.show("sticks","sele")
dist ////60/CA/,////61/CA/
select sele, resi 60+61
cmd.show("sticks","sele")
dist ////68/CA/,////71/CA/
select sele, resi 68+71
cmd.show("sticks","sele")
dist ////83/CA/,////84/CA/
select sele, resi 83+84
cmd.show("sticks","sele")
dist ////22/CA/,////118/CA/
select sele, resi 22+118
cmd.show("sticks","sele")
dist ////75/CA/,////78/CA/
select sele, resi 75+78
cmd.show("sticks","sele")
dist ////61/CA/,////63/CA/
select sele, resi 61+63
cmd.show("sticks","sele")
dist ////38/CA/,////41/CA/
select sele, resi 38+41
cmd.show("sticks","sele")
dist ////91/CA/,////103/CA/
select sele, resi 91+103
cmd.show("sticks","sele")
dist ////36/CA/,////37/CA/
select sele, resi 36+37
cmd.show("sticks","sele")
dist ////123/CA/,////129/CA/
select sele, resi 123+129
cmd.show("sticks","sele")
dist ////60/CA/,////62/CA/
select sele, resi 60+62
cmd.show("sticks","sele")
dist ////97/CA/,////113/CA/
select sele, resi 97+113
cmd.show("sticks","sele")
dist ////52/CA/,////53/CA/
select sele, resi 52+53
cmd.show("sticks","sele")
dist ////84/CA/,////90/CA/
select sele, resi 84+90
cmd.show("sticks","sele")
dist ////71/CA/,////73/CA/
select sele, resi 71+73
cmd.show("sticks","sele")
dist ////135/CA/,////137/CA/
select sele, resi 135+137
cmd.show("sticks","sele")
dist ////97/CA/,////101/CA/
select sele, resi 97+101
cmd.show("sticks","sele")
dist ////73/CA/,////74/CA/
select sele, resi 73+74
cmd.show("sticks","sele")
dist ////56/CA/,////94/CA/
select sele, resi 56+94
cmd.show("sticks","sele")
dist ////46/CA/,////90/CA/
select sele, resi 46+90
cmd.show("sticks","sele")
dist ////46/CA/,////50/CA/
select sele, resi 46+50
cmd.show("sticks","sele")
dist ////22/CA/,////119/CA/
select sele, resi 22+119
cmd.show("sticks","sele")
dist ////128/CA/,////130/CA/
select sele, resi 128+130
cmd.show("sticks","sele")
dist ////11/CA/,////107/CA/
select sele, resi 11+107
cmd.show("sticks","sele")
dist ////121/CA/,////131/CA/
select sele, resi 121+131
cmd.show("sticks","sele")
dist ////40/CA/,////41/CA/
select sele, resi 40+41
cmd.show("sticks","sele")
dist ////46/CA/,////83/CA/
select sele, resi 46+83
cmd.show("sticks","sele")
dist ////121/CA/,////129/CA/
select sele, resi 121+129
cmd.show("sticks","sele")
dist ////1/CA/,////28/CA/
select sele, resi 1+28
cmd.show("sticks","sele")
dist ////49/CA/,////56/CA/
select sele, resi 49+56
cmd.show("sticks","sele")
dist ////132/CA/,////135/CA/
select sele, resi 132+135
cmd.show("sticks","sele")
dist ////67/CA/,////86/CA/
select sele, resi 67+86
cmd.show("sticks","sele")
bond ////70/CA/,////109/CA/
bond ////79/CA/,////102/CA/
bond ////31/CA/,////66/CA/
select sele, resi 13
cmd.show("dots","sele")
cmd.show("sticks","sele")
set dot_density, 2
set dot_radius, 0
set dot_width, 1
cmd.hide("lines","ARSC_plmDCA")
set_bond line_color,white, ////70/CA/,////109/CA/
show lines, ////70/CA/ | ////109/CA/
set_bond line_color,white, ////79/CA/,////102/CA/
show lines, ////79/CA/ | ////102/CA/
set_bond line_color,white, ////31/CA/,////66/CA/
show lines, ////31/CA/ | ////66/CA/
cmd.show("ribbon","ARSC_plmDCA")
cmd.hide("((byres (ARSC_plmDCA))&(n. c,o,h|(n. n&!r. pro)))")
cmd.color(13,"dist*")
cmd.hide("labels","all")
cmd.bg_color('gray')
set ortho,0
view 1, store
set depth_cue = 0
set ray_trace_fog, 0
set ray_shadows, 0
ray 1000,1200
png ARSC_plmDCA.png, width=1000, height=1200, dpi=300, ray=1