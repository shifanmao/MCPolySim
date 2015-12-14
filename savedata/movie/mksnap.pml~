delete *
reset

cd pdb

#movie.load snap*.pdb,snap
load snap010.pdb,snap
rotate y,90,state=-1
rotate z,90,state=-1

#run ../color_b.py
#color_b snap,mode=hist,gradient=bgr,minimum=0,maximum=1,nbins=40,sat=1.,value=1.

hide all
#hide (name A2)
#hide (name A3)

show (name A1)
color blue, (name A1)
#set stick_radius=0.2
#show sticks, (name A1)
alter (name A1),vdw=0.5
show spheres, (name A1)

show (name A2)
color red, (name A2)
#set stick_radius=0.2
#show sticks, (name A2)
alter (name A2),vdw=0.5
show spheres, (name A2)


#show (name A3)
#color blue, (name A3)
#set stick_radius=0.2
#show sticks, (name A3)

bg_color white
reset

cd ..

#cd ../png
#set ray_trace_frames = 1
#mpng snap
