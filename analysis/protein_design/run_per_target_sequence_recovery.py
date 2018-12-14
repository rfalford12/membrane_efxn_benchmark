#!/usr/bin/env python
# @file: run_per_target_sequence_recovery.py
# @brief: Calculate sequence recovery over individual targets
# @author: Rebecca F. Alford (ralford3@jhu.edu)

# Main function that reads in a list of native targets and designed targets
# hyd values in each lipid composition should live in the inputs directory
# will also need to know which lipid compositions to compare

# A function that compares residues and atoms with hyd values
# If there is an atom that has a different hyd value, re-compare the residue
# (may need to round the hyd value to the nearest ??)

# Then, make a true/false vector as - should we measure recovery at this residue? 
# then then - use this script to measure sequence recovery and design recovery and then print to a file

import sys, os
import numpy as np
from pyrosetta import *

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )


workdir = "/Users/ralford/research/membrane_efxn_benchmark/analysis"
comp1_file = "1AFO_tr_DLPC_hydration.dat"
comp2_file = "1AFO_tr_DOPC_hydration.dat" 

def classify_hydration( h ): 

	if ( h < 0.25 ):
		return "lipid"
	elif ( h >= 0.25 and h < 0.75 ): 
		return "interfacial" 
	else: 
		return "aqueous"

def read_per_res_hydration( fn ): 

	# Read the data for amino acid composition 1
	with open( fn, 'rb') as f: 
		hydration_vals = f.readlines()
		hydration_vals = [ x.strip() for x in hydration_vals  ]
		hydration_vals  = [ x.split() for x in hydration_vals  ]

	per_res_hydration = []
	curr_resnum = hydration_vals[1][0]
	last_resnum = curr_resnum
	hatoms_hydration = []
	for i in range(1, len(hydration_vals)):
		if ( last_resnum == curr_resnum ): 
			hval = round( float( hydration_vals[i][2] ), 2 )
			hatoms_hydration.append( classify_hydration( hval ) )
			last_resnum = hydration_vals[i][0]
		else: 
			per_res_hydration.append( hatoms_hydration )
			hatoms_hydration = []
			hval = round( float( hydration_vals[i][2] ), 2 )
			hatoms_hydration.append( classify_hydration( hval ) )
			curr_resnum = hydration_vals[i][0]
			last_resnum = curr_resnum

	return per_res_hydration

def calculate_subset_with_different_hyd( fn_comp1, fn_comp2 ): 

	# Read hydration data for both lipid compositions
	composition1_hyd = read_per_res_hydration( fn_comp1 )
	composition2_hyd = read_per_res_hydration( fn_comp2 )
	assert len(composition1_hyd) == len(composition2_hyd)
	protein_length = len(composition1_hyd)

	subset_residues = []
	for res in range(0, protein_length): 

		assert len(composition1_hyd[res]) == len(composition2_hyd[res])
		res_length = len(composition1_hyd[res])

		for atom in range(1, res_length): 

			hyd_comp1 = composition1_hyd[res][atom]
			hyd_comp2 = composition2_hyd[res][atom]

			if ( hyd_comp1 != hyd_comp2 ): 
				subset_residues.append( res )

	# Convert subset residues into a nonredundant set
	subset_residues_set = set(subset_residues)
	return subset_residues_set

def compute_sequence_recovery_on_subset( native_pose, design_pose, subset ): 

	n_canonical = 20
	n_correct = np.zeros( n_canonical )
	n_native = np.zeros( n_canonical )
	n_designed = np.zeros( n_designed )

	# Log the native subset sequence
	native_sequence = []
	for r in range(0, native_pose.total_residue()): 
		if ( r in subset and native_pose.residue(r).is_protein() ): 
			native_sequence.append( native_pose.residue(r).name1() )

	# Calculate the n_native array
	for r in subset: 
		n_native[ native_pose.residue(r).aa() ]++; 


	design_sequence = []
	for r in range(0, design_pose.total_residue()): 

		# Measure sequence recovery
		n_designed[ design_pose.residue(r).aa() ]++
		if ( native_pose.redesign(r).aa() == design_pose.residue(r).aa() ): 
			n_correct[ design_pose.residue(r).aa() ]++


def main( argv ): 

	parser = OptionParser(usage="usage %prog --native_pdb_list natives.list --redesign_pdb_list redesign.list" )
	parser.set_description(main.__doc__)

	parser.add_option('--native_pdb_list', '-n', 
		action="store", 
		help="Name of file containing native PDBs",)

	parser.add_option('--redesign_pdb_list', '-r', 
		action="store", 
		help="Name of file containing redesign PDBs",)

	parser.add_option('--redesign_composition1', '-a'. 
		action="store", 
		help="Name of the file containing hydration value data files for c1", )

	parser.add_option('--redesign_composition2', '-b', 
		action="store", 
		help="Name of the file containing hydration value data files for c2", )

	parser.add_option('--composition1', '-l', 
		action="store", 
		help="Name of lipid compositoin 1", )

	parser.add_option('--composition2', '-m', 
		action="store", 
		help-"Name of lipid composition 2", )

	parser.add_option

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Read the native PDB list
    with open( Options.native_pdb_list, 'rb' ) as natives:
    	native_pdbs = natives.readlines()
    	native_pdbs = [ x.strip() for x in native_pdbs ]

    # Read the redesign PDB list
    with open( Options.redesign_pdb_list, 'rb' ) as redesigned:
    	redesgin_pdbs = redesgined.readlines()
    	redesgin_pdbs = [ x.strip() for x in redesign_pdbs ]

    # Assert the lists have the same length
    assert len(native_pdbs) == len(redesign_pdbs)
    n_targets = len(native_pdbs)

    # Read in the native poses
    native_poses = []
    for pdb in native_pdbs: 
    	pose = pose_from_pdb( pdb )
    	native_poses.append( pose )

    # Read in the redesigned poses from composition 1
    redesigned_poses = []
    for pdb in redesign_pdbs: 
    	pose = pose_from_pdb( pdb )
    	redesigned_poses.append( pose )

    # Read in the redesigned poses from composition 2

    # Read in the hydration parameters from lipid composition 1
    # Read in the hydration parameters from lipid composition 2

    # For each target, calculate recovery over the residue subsets
    for t in range(0, n_targets): 

    	# Find the mutual set of residues with a different hyd value
    	residue_subset = calculate_subset_with_different_hyd( target_composition1[t], target_composition2[t] )
		
    	# Calculate sequence recovery for poses designed in lipid composition 1
		n_correct_c1, n_native_c1, n_designed_c1 = compute_sequence_recovery_on_subset( native_poses[t], redesigned_in_c1_poses[t], residue_subset )

		# Calculate sequence recovery for poses designed in lipid composition 2
		n_correct_c2, n_native_c2, n_designed_c2 = compute_sequence_recovery_on_subset( native_poses[t], redesigned_in_c2_poses[t], residue_subset )


subset_residues_set = calculate_subset_with_different_hyd( comp1_file, comp2_file )




if __name__ == "__main__" : main(sys.argv)










