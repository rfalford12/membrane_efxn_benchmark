# @file: predict_ddG_of_insertion.py
# @brief: Calculate the ddG of inserting a peptide into the bilayer from solution
# @author: Rebecca Alford (ralford3@jhu.edu)

#from pyrosetta import *
from string import Template

import sys, os
import commands
import random

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main( args ):

    # Read options from the commandline
    parser = OptionParser( usage="usage: %prog --rscript landscape.xml --pdb test.pdb --span span.pdb --energy_fxn mpframework_fa_2007" )
    parser.set_description(main.__doc__)

    parser.add_option('--rosettaexe', '-r',
        action="store",
        help="Path to Rosetta executable", )

    parser.add_option('--rscript', '-t',
        action="store",
        help="Rosetta script used to perform grid search calculations",)

    parser.add_option('--pdb', '-i',
        action="store",
        help="PDB file containing coordinates for peptide",)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Energy function to use",)

    parser.add_option('--outfile', '-o',
        action="store",
        help="Output file containing ddg data"
    )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    if ( not Options.pdb or not Options.energy_fxn ):
        print "Missing required option --pdb or --energy_fxn! Exiting..."
        sys.exit()

    if ( not Options.rscript ):
        print "Missing required path to grid search Rosetta script"
        sys.exit()

    # Make a directory for case-specific data
    casestr = Options.pdb.split("/")[-1]
    case = casestr[:-4]
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + Options.energy_fxn + "/ddG-of-insertion/" + case
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Execute C++ task from this python script
    xml_script = Options.rscript
    s = Template( " -overwrite -parser:script_vars sfxn_weights=$sfxn_weights -parser:protocol $rscript -in:file:s $pdb -mp:setup:spanfiles single_TM_mode -out:path:all $outdir" )
    arguments = s.substitute( rscript=xml_script, pdb=Options.pdb, sfxn_weights=Options.energy_fxn, outdir=outdir )
    os.system( Options.rosettaexe + " " + arguments )

    # Check that the output file is present
    if ( not os.path.isfile( outdir + "/" +  case + "_" + Options.energy_fxn + "_landscape.dat" ) ):
        print "Energy landscpae calculations did not complete! Exiting..."
        sys.exit()

    # Read file into 3D data structure
    landscapefile = outdir + "/" + case + "_" + Options.energy_fxn + "_landscape.dat"
    with open( landscapefile, 'r' ) as f:
        content = f.readlines()
    content = [ x.strip() for x in content ]
    content = [ x.split(' ') for x in content ]

    # Read lines into an array of data triplets
    zcoords = []
    angles = []
    scores = []
    for x in content:
        if ( x[0] != "zcoord" ):
            zcoords.append( float(x[0]) )
            angles.append( float(x[1]) )
            scores.append( float(x[2]) )

    # Get the energy of the non-inserted pose
    max_index = zcoords.index( max(zcoords) )
    non_inserted_score = scores[ max_index ]

    # Get the lowest energy coordinate
    min_index = scores.index( min(scores) )
    lowest_zcoord = zcoords[ min_index ]
    lowest_angle = angles[ min_index ]
    lowest_score = scores[ min_index ]
    ddG_from_lowest = lowest_score - non_inserted_score

    # Get the score of the membrane centered pose ("fully inserted")
    rounded_zcoord = [ round(x) for x in zcoords ]
    inserted_index = rounded_zcoord.index( 0 )
    inserted_zcoord = zcoords[ inserted_index ]
    inserted_angle = angles[ inserted_index ]
    inserted_score = scores[ inserted_index ]
    ddG_from_inserted = inserted_score - non_inserted_score

    with open( Options.outfile, 'a' ) as f:
        s = Template( "$case $non_inserted_score $lowest_zcoord $lowest_angle $lowest_score $ddG_from_lowest $inserted_zcoord $inserted_angle $inserted_score $ddG_from_inserted\n" )
        outstr = s.substitute( case=case, non_inserted_score=non_inserted_score, lowest_zcoord=lowest_zcoord, lowest_angle=lowest_angle, lowest_score=lowest_score, ddG_from_lowest=ddG_from_lowest, inserted_zcoord=inserted_zcoord, inserted_angle=inserted_angle, inserted_score=inserted_score, ddG_from_inserted=ddG_from_inserted )
        f.write( outstr )

if __name__ == "__main__" : main(sys.argv)
