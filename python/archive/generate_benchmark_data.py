#!/usr/bin/env python
# @file: generate_benchmark_data.py
# @brief: Master script for running membrane energy function benchmarks on Jazz
# @notes: This script should **always** be run from the home membrane-efxn directory
# @author: Rebecca F. Alford (rfalford12@gmail.com)

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

##############################################################################
### Global data pertaining to the benchmark run
benchmark = "/home/ralford/membrane-efxn-benchmark/"
rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
rosettadir_stable = "/home/ralford/apps/Rosetta-stable/main/source/bin/"
buildenv = "release"
compiler = "gcc"
##############################################################################

def write_and_submit_condor_script( path, name, executable, arguments, queue_no=1 ):

    # Given a test path and
    filename = path + "/" + name + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + name + ".out\n" )
        f.write( "error = " + path + "/" + name + ".err\n" )
        f.write( "notify_user = rfalford12@gmail.com\n" )
        f.write( "request_memory = 15000\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

## Intiution Test: Compute the 2D membrane energy landscape of a monomeric alpha helix
def test_monomer_landscape( energy_fxn, rosetta_exe_path ):

    print "Initializing test: monomer-landscape"

    # Read list of monomer landscape test cases
    path_to_test = benchmark + "tests/intuition/test-monomer-landscape"
    list_of_test_cases = path_to_test + "/inputs/helices.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate path to executable
    executable = rosetta_exe_path + "rosetta_scripts" + ".linux" + compiler + buildenv
    xml_script = path_to_test + "/monomer-landscape.xml"

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/monomer-landscape"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = path_to_test + "/inputs/" + case + "/" + case + ".pdb"
        spanfile = path_to_test + "/inputs/" + case + "/" + case + ".span"
        arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting monomer-landscape test case:", case
        write_and_submit_condor_script( outdir, case, executable, arguments )

def test_aro_landscape( energy_fxn, rosetta_exe_path ):

    print "Initializing test: aro-landscape"

    # Read list of aro capped helix landscape test cases
    path_to_test = benchmark + "tests/intuition/test-aro-landscape"
    list_of_test_cases = path_to_test + "/inputs/aro_helices.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to executable
    executable = rosetta_exe_path + "rosetta_scripts.linux" + compiler + buildenv
    xml_script = path_to_test + "/aro-landscape.xml"

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/aro-landscape"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = path_to_test + "/inputs/" + case
        arguments = "-overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles from_structure -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script

        # Generate a condor submission file and submit the job to Jazz
        condor_case_name = case.split("/")[1]
        print "Submitting aro-landscape test case:", condor_case_name
        write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

def test_lk_landscape( energy_fxn, rosetta_exe_path ):

    print "Initializing test: lk-landscape"

    # Read list of aro capped helix landscape test cases
    path_to_test = benchmark + "tests/intuition/test-lk-landscape"
    list_of_test_cases = path_to_test + "/inputs/lk_peptides.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to executable
    executable = rosetta_exe_path + "rosetta_scripts.linux" + compiler + buildenv
    xml_script = path_to_test + "/lk-landscape.xml"

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/lk-landscape"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = path_to_test + "/inputs/" + case
        arguments = "-overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles from_structure -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script

        # Generate a condor submission file and submit the job to Jazz
        condor_case_name = case.split("/")[1]
        print "Submitting lk-landscape test case:", condor_case_name
        write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

def test_ddG_of_mutation( energy_fxn ):

    ## TODO: Going to update this to do a flexible backbone based calculation?? 

    print "Initializing test: ddG-of-mutation"

    # Specify path to test, setup new directories
    path_to_test = benchmark + "tests/validation/test-ddG-of-mutation"
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/ddG-of-mutation"
    os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    # Specify path to python script and arguments
    python_script = path_to_test + "/predict_ddG.py"
    energy_function = energy_fxn
    mlist = path_to_test + "/inputs/all_mutations.dat"
    s = Template( "--energy_fxn $energy_func --mutation_list $list_of_mutations --outdir  $outdir")
    arguments = s.substitute( energy_func=energy_function, list_of_mutations=mlist, outdir=outdir )

    # Execute the command within the virtual environment
    print "Submitting all ddG of mutation test cases"
    os.system( "python " + python_script + " " + arguments )

def test_ddG_of_insertion( energy_fxn, rosetta_exe_path ):

    print "Initializing test: ddG-of-insertion"

    # Read list of aro capped helix landscape test cases
    path_to_test = benchmark + "test/validation/test-ddG-of-insertion"
    list_of_test_cases = path_to_test + "/inputs/"

    # Read list of aro capped helix landscape test cases
    path_to_test = benchmark + "tests/validation/test-ddG-of-insertion"
    list_of_test_cases = path_to_test + "/inputs/insertion_peptide.dat"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to executable
    executable = rosetta_exe_path + "rosetta_scripts.linux" + compiler + buildenv
    xml_script = path_to_test + "/landscape.xml"

    # Path to test script
    test_script = path_to_test + "/predict_ddG_of_insertion.py"

    # Change directories to a data analysis dir
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/ddG-of-insertion"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Make a master file in the output directory that logs positions
    filename = outdir + "/ddG-of-insertion-results.dat"
    with open( filename, 'a' ) as f:
        f.write( "case non_inserted_score lowest_zcoord lowest_angle lowest_score ddG_from_lowest inserted_zcoord inserted_angle inserted_score ddG_from_inserted\n" )

        # For each test case, generate specific arguments, condor_files, and then run
    for case in test_cases:

        case = case.split(' ')
        if ( case[0] != "name" ):

            # Setup case-specific variables (pdbfile, spanfile, xmlargs)
            pdbfile = path_to_test + "/inputs/" + case[0] + "/" + case[0] + "_startstruc.pdb"
            s = Template( " --pdb $pdb --rscript $xml --energy_fxn $sfxn --rosettaexe $rosetta --outfile $outfile" )
            arguments = s.substitute( pdb=pdbfile, xml=xml_script, sfxn=Options.energy_fxn, rosetta=executable, outfile=filename )

            # Run Python script
            print "python " + test_script + arguments
            os.system( "python " + test_script + arguments )

def test_pH_dependent_insertion( energy_fxn, rosetta_exe_path ):

    print "Initializing test: pH-dependent-insertion"

    # Read in library of pH sensitive insertion peptides
    path_to_test = benchmark + "tests/validation/test-pH-dependent-insertion"
    list_of_test_cases = path_to_test + "/inputs/pH-inserted-peptides.dat"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/pH-dependent-insertion"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Generate a path to executable
    executable = rosetta_exe_path + "rosetta_scripts.linux" + compiler + buildenv
    xml_script = path_to_test + "/pH-dependent-landscape.xml"

    # Path to test script
    test_script = path_to_test + "/predict_ddG_of_pH_sensitive_insertion.py"

    # Make a master file in the output directory that logs positions
    filename = outdir + "/pH-sensitive-insertion-results.dat"
    with open( filename, 'a' ) as f:
        f.write( "case non_inserted_score lowest_zcoord lowest_angle lowest_score ddG_from_lowest inserted_zcoord inserted_angle inserted_score ddG_from_inserted\n" )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test_cases:

        case = case.split(' ')
        if ( case[0] != "Sequence" ):

            # Setup case-specific variables (pdbfile, spanfile, xmlargs)
            pdbfile = path_to_test + "/inputs/" + case[0] + ".pdb"
            s = Template( " --pdb $pdb --rscript $xml --energy_fxn $sfxn --rosettaexe $rosetta --outfile $outfile" )
            arguments = s.substitute( pdb=pdbfile, xml=xml_script, sfxn=Options.energy_fxn, rosetta=executable, outfile=filename )

            # # Run Python script
            os.system( "python " + test_script + arguments )

def test_tilt_angle( energy_fxn, rosetta_exe_path ):

    print "Initializing test: tilt-angle"

    # This test is a head-fake: it uses the same dataset as the monomer landscape test
    # Thus, we just copy the data over into a separate file

    # Change directories to a data analysis dir
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/tilt-angle"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Copy data into a separate folder for analysis
    origin = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/monomer-landscape"
    copy_cmd = "cp -r " + origin + " ."
    os.system( copy_cmd )

def test_sequence_recovery( energy_fxn, rosetta_exe_path ):

    print "Initializing test: sequence-recovery"

    # Read list of monomeric TM spanning chains
    path_to_test = benchmark + "tests/validation/test-seq-recovery"
    list_of_test_cases = path_to_test + "/inputs/monomer_chains.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to the fixbb executable
    executable = rosetta_exe_path + "fixbb.linux" + compiler + buildenv

    # Change directories into a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/seq-recovery"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/seq-recovery/" + case
        os.system( "mkdir " + outdir )
        os.chdir( outdir )

        # Setup arguments by substitution
        pdbfile = path_to_test + "/inputs/" + case + "/" + case + "_tr_ignorechain.pdb"
        spanfile = path_to_test + "/inputs/" + case + "/" + case + "_tr.span"
        s = Template( "-in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res -restore_talaris_behavior" )
        arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting sequence recovery test case:", case
        queue_no = 1
        write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )

def test_docking( energy_fxn, rosetta_exe_path ):

    print "Initializing test: docking"

    # ### Setup some general variable
    path_to_test = benchmark + "tests/validation/test-docking"
    executable = rosetta_exe_path + "mp_dock.linux" + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/docking"
    os.system( "mkdir " + base_outdir )

    ### Test #1: Docking set from Lomize et al. 2017 (small homodimers)
    list_of_test01_cases = path_to_test + "/inputs/small-homodimer-set/dimers.list"
    with open( list_of_test01_cases, 'rb' ) as f:
        test01_cases = f.readlines()
    test01_cases = [ x.strip() for x in test01_cases ]

    test01_outdir = benchmark + "data/" + energy_fxn + "/docking/small-homodimer-set"
    os.system( "mkdir " + test01_outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test01_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/docking/small-homodimer-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup arguments by substitution
        native = path_to_test + "/inputs/small-homodimer-set/" + case + "/" + case + "_AB_tr.pdb"
        prepacked = path_to_test + "/inputs/small-homodimer-set/" + case + "/" + case + "_AB_tr.prepack.pdb"
        spanfile = path_to_test + "/inputs/small-homodimer-set/" + case + "/" + case + "_AB.span"
        s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners A_B -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir" )
        arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting docking test case from small homodimer set:", case
        queue_no = 300
        write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )

    ### Test #2: Docking set from Alford & Koehler Leman et al. 2015 (large homodimers)
    list_of_test02_cases = path_to_test + "/inputs/large-homodimer-set/dimers.list"
    with open( list_of_test02_cases, 'rb' ) as f:
        test02_cases = f.readlines()
    test02_cases = [ x.strip() for x in test02_cases ]

    test02_outdir = benchmark + "data/" + energy_fxn + "/docking/large-homodimer-set"
    os.system( "mkdir " + test02_outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test02_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/docking/small-homodimer-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup arguments by substitution
        native = path_to_test + "/inputs/large-homodimer-set/" + case + "/" + case + "_AB_tr.pdb"
        prepacked = path_to_test + "/inputs/large-homodimer-set/" + case + "/" + case + "_AB_tr.prepack.pdb"
        spanfile = path_to_test + "/inputs/large-homodimer-set/" + case + "/" + case + "_AB.span"
        s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners A_B -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir" )
        arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting docking test case from large homodimer set:", case
        queue_no = 300
        write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )


def test_decoy_discrimination( energy_fxn, rosetta_exe_path ):

    print "Initializing test: decoy-discrimination"

    ### Setup some general variables
    path_to_test = benchmark + "tests/validation/test-decoy-discrimination"
    executable = rosetta_exe_path + "rosetta_scripts.linux" + compiler + buildenv
    xml_script = path_to_test + "/mp_refinement.xml"
    base_outdir = benchmark + "data/" + energy_fxn + "/decoy-discrimination"
    os.system( "mkdir " + base_outdir)

    ### Test Set #1: Yarov-Yaravoy Low Resolution Decoys (Membrane ab initio generated)
    list_of_test01_cases = path_to_test + "/inputs/yarov-yaravoy-set/decoy_sets.list"
    with open( list_of_test01_cases, 'rb' ) as f:
        test01_cases = f.readlines()
    test01_cases = [ x.strip() for x in test01_cases ]

    outdir_test01 = benchmark + "data/" + energy_fxn + "/decoy-discrimination/yarov-yaravoy-set"
    os.system( "mkdir " + outdir_test01 )
    os.chdir( outdir_test01 )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test01_cases:

        outdir = benchmark + "data/" + energy_fxn + "/decoy-discrimination/yarov-yaravoy-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        for i in range(1, 51):

            # Setup case-specific variables (pdbfile, spanfile, xmlargs)
            native = path_to_test + "/inputs/yarov-yaravoy-set/" + case + "/" + case + "_native.pdb"
            spanfile = path_to_test + "/inputs/yarov-yaravoy-set/" + case + "/" + case + ".span"
            modelslist = path_to_test + "/inputs/yarov-yaravoy-set/" + case + "/models." + str(i) + ".list"
            s = Template( "-cddockingrestore_talaris_behavior -overwrite -in:file:native $native -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
            arguments = s.substitute( modellist=modelslist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)

            # Generate a condor submission file and submit the job to Jazz
            condor_case_name = case + "_models_" + str(i)
            print "Submitting decoy-discrimination test case from Yarov-Yaravoy set:", condor_case_name
            write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

    ### Test Set #2: Dutagaci High Resolution Decoys (Molecular dynamics generated)
    list_of_test02_cases = path_to_test + "/inputs/dutagaci-set/decoy_sets.list"
    with open( list_of_test02_cases, 'rb' ) as f:
        test02_cases = f.readlines()
    test02_cases = [ x.strip() for x in test02_cases ]

    outdir_test02 = benchmark + "data/" + energy_fxn + "/decoy-discrimination/dutagaci-set"
    os.system( "mkdir " + outdir_test02 )
    os.chdir( outdir_test02 )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test02_cases:

        outdir = benchmark + "data/" + energy_fxn + "/decoy-discrimination/dutagaci-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        native = path_to_test + "/inputs/dutagaci-set/" + case + "/" + case + "_native.pdb"
        spanfile = path_to_test + "/inputs/dutagaci-set/" + case + "/" + case + ".span"
        modelslist = path_to_test + "/inputs/dutagaci-set/" + case + "/decoys.list"
        s = Template( "-overwrite -in:file:native $native -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
        arguments = s.substitute( modellist=modelslist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)

        # Generate a condor submission file and submit the job to Jazz
        condor_case_name = case + "_dutagaci"
        print "Submitting decoy-discrimination test case from Dutagaci set:", condor_case_name
        write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_t01 --which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option('--which_tests', '-w',
        action="store",
        help="Which tests should I run? all, intiution, or validation", )

    parser.add_option('--stable_master', '-s', 
        action="store", 
        help="For reference runs, use the stable Rosetta master branch",
        )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    if ( not Options.energy_fxn or not Options.which_tests ):
        print "Missing required options --energy_fxn and --which_tests! Exiting..."
        sys.exit()

    if ( not Options.stable_master ): 
        print "Missing required option --stable_master"
        sys.exit()

    if ( Options.which_tests != 'all' and Options.which_tests != 'intuition' and Options.which_tests != 'validation' ):
        print "Illegal value for option --which_tests! Found", Options.which_tests
        sys.exit()

    # Choose rosetta path based on input options
    rosetta_exe_path = ""
    if ( Options.stable_master == "true" ): 
        rosetta_exe_path = rosettadir_stable
    else: 
        rosetta_exe_path = rosettadir

    # Make an output data directory
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + Options.energy_fxn
    os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    if ( Options.which_tests == 'all' or Options.which_tests == 'intuition' ):
        print "Initializing intiution Tests"
        #test_monomer_landscape( Options.energy_fxn, rosetta_exe_path )
        #test_aro_landscape( Options.energy_fxn, rosetta_exe_path )
        #test_lk_landscape( Options.energy_fxn, rosetta_exe_path )

    if ( Options.which_tests == 'all' or Options.which_tests == 'validation' ):
        print "Initializing validation tests"
        #test_docking( Options.energy_fxn, rosetta_exe_path )
        test_sequence_recovery( Options.energy_fxn, rosetta_exe_path )
        #test_decoy_discrimination( Options.energy_fxn, rosetta_exe_path )
        #test_ddG_of_mutation( Options.energy_fxn )
        #test_ddG_of_insertion( Options.energy_fxn, rosetta_exe_path )
        #test_pH_dependent_insertion( Options.energy_fxn, rosetta_exe_path )
        #test_tilt_angle( Options.energy_fxn, rosetta_exe_path )


if __name__ == "__main__" : main(sys.argv)
