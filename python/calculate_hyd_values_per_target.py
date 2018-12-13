#!/usr/bin/env python
# @file: calculate_hyd_values_per_target.py
# @brief: Calculate hydration values per different targets
# @author: Rebecca F. Alford (ralford3@jhu.edu)

##############################################################################
### Global data for benchmark runs on Jazz
rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
platform = "linux"
buildenv = "release"
compiler = "gcc"
##############################################################################


import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main( argv ): 

	parser = OptionParser(usage="usage %prog --pdb_list redesigned.list" )
	parser.set_description(main.__doc__)

	parser.add_option('--pdb_list', '-p', 
		action="store", 
		help="Name of file containing list of redesigned pdbs",)

	parser.add_option('--lipid_type', '-l', 
		action="store", 
		help="Name of lipid type",)

   	(options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    with open( Options.pdb_list, 'rb' ) as f: 
    	targets = f.readlines()
    	targets = [ x.strip() for x in targets ]

    for t in range(0, len(targets)): 
		executable = rosetta_exe_path + "compute_lipid_specific_hyd" + "." + platform + compiler + buildenv
    	pdbfile = targets[t]
        spanfile = targets[t].split('_')[0] + "_tr.span"
        s = Template( " -overwrite -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -read_only_ATOM_entries -mp:lipids:temperature 37 -mp:lipids:composition $comp" )
        arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, comp=Options.lipid_type )
        os.system( executable + arguments )


if __name__ == "__main__" : main(sys.argv)


