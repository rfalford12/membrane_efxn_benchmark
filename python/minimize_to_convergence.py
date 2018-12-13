#!/usr/bin/env python
#@file: minimize_to_convergence.py
#@brief: Minimize ZMPSTE24 until convergence
#@author: Rebecca Alford (ralford3@jhu.edu)

from pyrosetta import *
from pyrosetta.teaching import *
init( extra_options=" -mp:lipids:use_implicit_lipids true -mp:lipids:composition DLPC -mp:lipids:temperature 20" )

# Read in PDB
pose = pose_from_file( "3gp6_A.pdb" )
add_memb = pyrosetta.rosetta.protocols.membrane.AddMembraneMover( "3gp6_A.span" )
add_memb.apply( pose )
sfxn = create_score_function( "franklin2018" ) #get_score_function()

# Going to do a quick side chain packing first
task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
pack_mover = PackRotamersMover(sfxn, task_pack)
pack_mover.apply(pose)
print("After pack score:" + str(sfxn.score(pose)))

min_mover = MinMover()
mm = MoveMap()
mm.set_bb(True)
mm.set_chi(True)
min_mover.movemap(mm)
min_mover.score_function(sfxn)

for i in range(0, 500):

    print("===================Round " + str(i) + "===================")
    min_mover.apply(pose)
    print (sfxn.score(pose))

pose.dump_pdb( "1qd6_min.pdb" )
