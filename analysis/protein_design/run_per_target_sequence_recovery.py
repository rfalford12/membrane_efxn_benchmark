#!/usr/bin/env python
# @file: run_per_target_sequence_recovery.py
# @brief: Calculate sequence recovery over individual targets
# @author: Rebecca F. Alford (ralford3@jhu.edu)

<<<<<<< HEAD
import sys, os
=======
import sys
import os
>>>>>>> e39cb737e9de9d1dd7d49aaedcecadfb3e87306c
import numpy as np
from pyrosetta import *

from optparse import OptionParser, IndentedHelpFormatter
<<<<<<< HEAD
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

workdir = "/Users/ralford/research/membrane_efxn_benchmark/analysis"
=======
_script_path_ = os.path.dirname(os.path.realpath(__file__))

workdir = "/Users/ralford/research/membrane_efxn_benchmark/analysis"
benchmark_data = "/home/ralford/membrane_efxn_benchmark/data/franklin2018/"
aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
>>>>>>> e39cb737e9de9d1dd7d49aaedcecadfb3e87306c

def classify_hydration(h):

	if (h < 0.25):
		return "lipid"
	elif (h >= 0.25 and h < 0.75):
		return "interfacial"
	else:
		return "aqueous"


def read_per_res_hydration(fn):

	# Read the data for amino acid composition 1
	with open(fn, 'rb') as f:
		hydration_vals = f.readlines()
		hydration_vals = [x.strip() for x in hydration_vals]
		hydration_vals = [x.split() for x in hydration_vals]

	per_res_hydration = []
	curr_resnum = hydration_vals[1][0]
	last_resnum = curr_resnum
	hatoms_hydration = []
	for i in range(1, len(hydration_vals)):
		if (last_resnum == curr_resnum):
			hval = round(float(hydration_vals[i][2]), 2)
			hatoms_hydration.append(classify_hydration(hval))
			last_resnum = hydration_vals[i][0]
		else:
			per_res_hydration.append(hatoms_hydration)
			hatoms_hydration = []
			hval = round(float(hydration_vals[i][2]), 2)
			hatoms_hydration.append(classify_hydration(hval))
			curr_resnum = hydration_vals[i][0]
			last_resnum = curr_resnum

	return per_res_hydration


def calculate_subset_with_different_hyd(fn_comp1, fn_comp2):

	# Read hydration data for both lipid compositions
	composition1_hyd = read_per_res_hydration(fn_comp1)
	composition2_hyd = read_per_res_hydration(fn_comp2)
	assert len(composition1_hyd) == len(composition2_hyd)
	protein_length = len(composition1_hyd)

	subset_residues = []
	for res in range(0, protein_length):

		c1_set = set(composition1_hyd[res])
		c2_set = set(composition2_hyd[res])
		if ( c1_set != c2_set ): 
			subset_residues.append(res+1)

	# Convert subset residues into a nonredundant set
	subset_residues_set = set(subset_residues)
	return subset_residues_set

def compute_sequence_recovery_on_subset( native_pose, design_pose, subset ):

	n_canonical = len(aas)
	n_correct = dict.fromkeys(aas)
	n_native = dict.fromkeys(aas)
	for i in n_correct: 
		n_correct[i] = 0.0
		n_native[i] = 0.0

	# Log the native subset sequence
	native_sequence = []
	for r in range(1, native_pose.total_residue()):
		if ( r in subset and native_pose.residue(r).is_protein() ):
			native_sequence.append( native_pose.residue(r).name1() )

	# Calculate the n_native array
	for r in subset:
	#for r in range(1, native_pose.total_residue()): 
		n_native[ native_pose.residue(r).name1() ] = n_native[ native_pose.residue(r).name1() ] + 1
		if ( native_pose.residue(r).name1() == design_pose.residue(r).name1() ): 
			n_correct[ design_pose.residue(r).name1() ] = n_correct[ design_pose.residue(r).name1() ] + 1

	# Calculate per-residue recovery
	recovery = {}
	total_native = 0
	total_correct = 0
	for aa in aas: 
		native = n_native[aa]
		correct = n_correct[aa]
		total_native += n_native[aa]
		total_correct += n_correct[aa] 

		if ( native == 0 ): 
			recovery[aa] = -1
		else: 
			recovery[aa] = round( correct/native, 3)

	total_recovered = round( total_correct/total_native, 3)

	return recovery, total_recovered, total_correct, total_native

def main(args):

	parser = OptionParser(
	    usage="usage %prog --native_pdb_list natives.list --redesign_pdb_list redesign.list")
	parser.set_description(main.__doc__)

	parser.add_option('--native_pdb_list', '-n',
		action="store",
		help="Name of file containing native PDBs",)

	parser.add_option('--redesign_pdb_list', '-a',
		action="store",
		help="Name of file containing redesigned PDBs",)

<<<<<<< HEAD
	parser.add_option('--composition1', '-l', 
=======
	parser.add_option('--composition1', '-l',
		action="store",
		help="Name of lipid compositoin 1", )

	parser.add_option('--composition2', '-m',
		action="store",
		help="Name of lipid composition 2", )

	parser.add_option('--output', '-o',
>>>>>>> e39cb737e9de9d1dd7d49aaedcecadfb3e87306c
		action="store", 
		help="Name of output filename", )

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options

	# Initialize PyRosetta
	init()

	# Read the different lipid composiitons
	lipid_type1 = Options.composition1
	lipid_type2 = Options.composition2

    # Read the native PDB list
	with open(Options.native_pdb_list, 'rt') as natives:
		native_pdbs = natives.readlines()
		native_pdbs = [x.strip() for x in native_pdbs]

<<<<<<< HEAD
    # Read the redesign PDB list for composition 1
    with open( Options.redesign_pdb_list, 'rb' ) as redesigned:
    	redesgin_pdbs = redesgined.readlines()
    	redesgin_pdbs = [ x.strip() for x in redesign_pdbs ]
=======
    # Read the redesign PDB list from lipid composition 1
	with open( Options.redesign_pdb_list, 'rt' ) as redesigned:
		designed_pdbs = redesigned.readlines()
		designed_pdbs = [ x.strip() for x in designed_pdbs ]
		redesign_pdbs_c1 = [ (benchmark_data + lipid_type1 + "/" + x) for x in designed_pdbs ]
		redesign_pdbs_c2 = [ (benchmark_data + lipid_type2 + "/" + x) for x in designed_pdbs ]
>>>>>>> e39cb737e9de9d1dd7d49aaedcecadfb3e87306c

    # Read hydration parameter files 

    # Assert the lists have the same length
	assert len(native_pdbs) == len(designed_pdbs)
	n_targets = len(native_pdbs)

    # Read in the native poses
	native_poses = []
	for pdb in native_pdbs: 
		pose = pose_from_pdb( pdb )
		native_poses.append( pose )

    # Read in the redesigned poses from composition 1
	redesigned_poses_c1 = []
	for pdb in redesign_pdbs_c1: 
		pose = pose_from_pdb( pdb )
		redesigned_poses_c1.append( pose )

 	# Read in the redesigned poses from composition 2
	redesigned_poses_c2 = []
	for pdb in redesign_pdbs_c2: 
		pose = pose_from_pdb( pdb )
		redesigned_poses_c2.append( pose )

	# Assign an output filename
	output_filename = Options.output
	with open( output_filename, 'wt' ) as f: 
		f.write( "target lipid_type aa recovery\n" )

	# For each target, calculate the difference in hydration, then the difference in recovery
	n_correct_c1 = 0
	n_native_c1 = 0
	n_correct_c2 = 0
	n_native_c2 = 0
	for i in range(0, n_targets): 

		# Get the target name
		target_name = [ x.split("_") for x in native_poses[i].pdb_info().name().split("/") ][6][0]

		# Read in the hydration parameters from lipid composition 1
		lipid_type_name1 = lipid_type1.split("_")[1]
		hydration_fn_c1 = benchmark_data + lipid_type1 + "/" + target_name + "_" + lipid_type_name1 + "_hydration.dat"
		lipid_type_name2 = lipid_type2.split("_")[1]
		hydration_fn_c2 = benchmark_data + lipid_type2 + "/" + target_name + "_" + lipid_type_name2 + "_hydration.dat"

		# Find the subset of residues that are affected by a change in hyd
		residue_subset = calculate_subset_with_different_hyd( hydration_fn_c1, hydration_fn_c2 )

		# Calculate the sequence recovery for proteins designed in lipid composition 1
		c1_recovery, c1_total_recovery, c1_correct, c1_native = compute_sequence_recovery_on_subset( native_poses[i], redesigned_poses_c1[i], residue_subset )
		n_correct_c1 += c1_correct
		n_native_c1 += c1_native

		# Calculate sequence recovery for poses designed in lipid composition 2
		c2_recovery, c2_total_recovery, c2_correct, c2_native = compute_sequence_recovery_on_subset( native_poses[i], redesigned_poses_c2[i], residue_subset )
		n_correct_c2 += c2_correct
		n_native_c2 += c2_native

		# Output the data for c1
		with open( output_filename, 'at' ) as f:
			# Recovery information for lipid composition 1
			f.write( target_name + " " + lipid_type_name1 + " total " + str(c1_total_recovery) + "\n"  )
			for aa in c1_recovery: 
				f.write( target_name + " " + lipid_type_name1 + " " + aa + " " + str(c1_recovery[aa]) + "\n" )

			# Recovery information for lipid composition 2
			f.write( target_name + " " + lipid_type_name2 + " total " + str(c2_total_recovery) + "\n"  )
			for aa in c2_recovery: 
				f.write( target_name + " " + lipid_type_name2 + " " + aa + " " + str(c2_recovery[aa]) + "\n" )

	overall_c1_recov = round( n_correct_c1 / n_native_c1, 3 )
	overall_c2_recov = round( n_correct_c2 / n_native_c2, 3 )
	print("Recovery in", lipid_type1, ":", overall_c1_recov, "with available positions n=", n_native_c1)
	print("Recovery in", lipid_type2, ":", overall_c2_recov, "with available positions n=", n_native_c2)

if __name__ == "__main__" : main(sys.argv)










