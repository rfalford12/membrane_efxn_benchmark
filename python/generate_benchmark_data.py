#!/usr/bin/env python
# @file: generate_benchmark_data_marcc.py
# @brief: Master script for running membrane energy function benchmarks on MARCC
# @notes: This script should **always** be run from the home membrane-efxn directory
# @author: Rebecca F. Alford (ralford3@jhu.edu)

## TODO: PH output needs to be put into separate folders otherwise the naming is ambiguous
## TODO: Change the paths to be compatible with this

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

##############################################################################
### Global data for benchmark runs on Jazz
#benchmark = "/home/ralford/membrane-efxn-benchmark/"
#rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
#rosettadir_stable = "/home/ralford/apps/Rosetta-stable/main/source/bin/"
#platform = "linux"
#buildenv = "release"
#compiler = "gcc"

##############################################################################
### Global data for benchmark runs on MARCC
benchmark = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/membrane_efxn_benchmark/"
rosettadir = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/Rosetta/main/source/bin/"
rosettadir_stable = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/Rosetta-stable/main/source/bin/"
platform = "mpi.linux" 
buildenv = "release"
compiler = "gcc"
##############################################################################


def write_and_submit_condor_script( path, name, executable, arguments, queue_no=1, high_mem=False ):

    # Given a test path and
    filename = path + "/" + name + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + name + ".out\n" )
        f.write( "error = " + path + "/" + name + ".err\n" )
        if ( high_mem == True ): 
            f.write( "request_memory = 15000\n")
        f.write( "notify_user = rfalford12@gmail.com\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

def write_and_submit_slurm_batch_script( path, name, executable, arguments, num_nodes=1 ): 

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + name + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on MARCC with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + name + "\n" )
        f.write( "#SBATCH --partition=parallel\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=60:0:0\n" )
        f.write( "#SBATCH --mem=120GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + name + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + name + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rfalford12@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        f.write( "module unload openmpi\n" )
        f.write( "module load intel-mpi\n" )

        # Provide a description of the job
        f.write( "ROSETTAEXE=" + executable + "\n" )
        f.write(  "echo Starting MPI job running $ROSETTAEXE\n" )

        # Run the job
        f.write( "mpirun $EXE " + arguments + "\n" )

        f.close()

    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    os.system( sbatch_command )

def run_energy_landscape_calc( energy_fxn, rosetta_exe_path, cluster_type, test_name, input_list, xml_protocol, restore, single_TM="false", pH="0" ): 
    """
    A general functions for running energy landscape calculations given a list of input helices

    Arguments:
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        test_name = Name of energy landscpae test variant
        input_list = List of input helices
        xml_protocol = Path to RosettaScript defining the landscape search protocol
        single_TM = use the single_TM_peptide spanfile option
    """

    print "Initializing energy landscape test for " + test_name

    # Read list of energy landscape test cases
    path_to_test = benchmark + "inputs/" + test_name 
    list_of_test_cases = path_to_test + "/" + input_list 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate path to executable
    executable = rosetta_exe_path + "rosetta_scripts" + "." + platform + compiler + buildenv
    xml_script = path_to_test + "/" + xml_protocol

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/" + test_name
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = path_to_test + "/inputs/" + case + "/" + case + ".pdb"

        # Should I use a single_TM_peptide span estimation or a user-provided spanfile? 
        spanfile = ""
        if ( single_TM == "false" ): 
            spanfile = path_to_test + "/inputs/" + case + "/" + case + ".span"
        else: 
            spanfile = "single_TM_mode"

        # Should I tune the pH of my simulation?
        arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script
        if ( pH != "0" ): 
            arguments = arguments + " -pH_mode true -value_pH " + pH
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting test case for " + test_name + ": " +  case
        if ( cluster_type == "SLURM" ): 
            write_and_submit_slurm_batch_script( outdir, case, executable, arguments, )
        else: 
            write_and_submit_condor_script( outdir, case, executable, arguments )

def run_ddG_of_mutation_test( energy_fxn, list_of_ddGs, test_name, restore ): 
    """
    A general function for calculating the ddG of single point mutations

    Arguments:  
        energy_fxn = energy function to run the protocol with
        list_of_ddGs = Relative path to list of ddGs to calculate
        test_name = Name of benchmark set for calculating ddGs of mutation
    """

    print "Initializing ddG-of-mutation test for set " + test_name

    # Specify path to test, setup new directories
    path_to_test = benchmark + "inputs/test_2.1_ddG_of_mutation"
    outdir = benchmark + "data/" + energy_fxn + "/test_2.1_ddG_of_mutation"
    os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    # Specify path to python script and arguments
    python_script = benchmark + "/python/test_2.1_predict_ddG.py"
    energy_function = energy_fxn
    mlist = path_to_test + "/" + list_of_ddGs
    s = Template( "--energy_fxn $energy_func --mutation_list $list_of_mutations --outdir  $outdir")
    arguments = s.substitute( energy_func=energy_function, list_of_mutations=mlist, outdir=outdir )
    if ( restore == True ): 
        arguments = arguments + " --restore"

    print "Submitting ddG of mutation test case " + test_name
    os.system( "python " + python_script + " " + arguments )

def run_fixed_backbone_design_calc( energy_fxn, rosetta_exe_path, cluster_type, restore ): 
    """
    A function for running the fixed backbone design calculations needed for the sequence recovery test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        restore = Restore to talalaris behavior prior to ref2015 for reference benchmark run
    """

    print "Initializing fixed backbone design calculations for sequence recovery test" 

    # Read list of test cases - typically full length transmembrane proteins
    inputs = benchmark + "inputs/test_3.1_sequence_recovery"
    list_of_test_cases = inputs + "/monomer_chains.list" 
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to the fixbb executable
    executable = rosetta_exe_path + "fixbb." + platform + compiler + buildenv

    # Change directories into a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/test_3.1_sequence_recovery"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.1_sequence_recovery/" + case
        os.system( "mkdir " + outdir )
        os.chdir( outdir )

        # Setup arguments by substitution
        pdbfile = inputs + "/" + case + "/" + case + "_tr_ignorechain.pdb"
        spanfile = inputs + "/" + case + "/" + case + "_tr.span"
        s = Template( " -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
        arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"

        # Generate a condor submission file and submit the job to Jazz
        print "Submitting fixed backbone design calculation for sequence recovery case:", case
        if ( cluster_type == "SLURM" ): 
            write_and_submit_slurm_batch_script( outdir, case, executable, arguments )
        else: 
            queue_no = 1
            high_mem = True 
            write_and_submit_condor_script( outdir, case, executable, arguments, queue_no, high_mem )

def run_docking_calc( energy_fxn, rosetta_exe_path, cluster_type, test_set, restore ): 
    """
    A function for running the docking calculations needed for the docking test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        test_set = Name of published set of test cases
        restore = Restore behavior to pre ref2015 for reference benchmark runs
    """

    print "Initializing test: docking"

    # Setup path to test sets and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "mp_dock" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup arguments by substitution
        native = inputs + "/" + test_set + "/" + case + "/" + case + "_AB_tr.pdb"
        prepacked = inputs + "/" + test_set + "/" + case + "/" + case + "_AB_tr.prepack.pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + "_AB.span"
        s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners A_B -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir" )
        arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"

        # Generate job submission file and then submit to cluster
        print "Submitting docking test case from set " + test_set + ":", case
        if ( cluster_type == "SLURM" ): 
            write_and_submit_slurm_batch_script( outdir, case, executable, arguments, str(10) )
        else: 
            queue_no = 300
            write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )

def run_decoy_discrimination_calc( energy_fxn, rosetta_exe_path, cluster_type, restore ): 
    """
    A function for running the docking calculations needed for the docking test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        restore = Option to restore talaris behavior before ref2015 for reference runs
    """

    print "Initializing test: decoy-discrimination"

    ### Setup some general variables
    inputs = benchmark + "inputs/test_3.3_decoy_discrimination"
    executable = rosetta_exe_path + "rosetta_scripts." + platform + compiler + buildenv
    xml_script = benchmark + "xml/test_3.3_decoy_refinement.xml"

    ### Make the base decoy discrimination output directory
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination"
    os.system( "mkdir " + base_outdir )

    ### Test Set #1: Yarov-Yaravoy Low Resolution Decoys (Membrane ab initio generated)
    list_of_test01_cases = inputs + "/yarov-yaravoy-set/decoy_sets.list"
    with open( list_of_test01_cases, 'rb' ) as f:
        test01_cases = f.readlines()
    test01_cases = [ x.strip() for x in test01_cases ]

    outdir_test01 = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/yarov-yaravoy-set"
    os.system( "mkdir " + outdir_test01 )
    os.chdir( outdir_test01 )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test01_cases:

        outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/yarov-yaravoy-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        for i in range(1, 51):

            # Setup case-specific variables (pdbfile, spanfile, xmlargs)
            native = inputs + "/yarov-yaravoy-set/" + case + "/" + case + "_native.pdb"
            spanfile = inputs + "/yarov-yaravoy-set/" + case + "/" + case + ".span"
            modelslist = inputs + "/yarov-yaravoy-set/" + case + "/models." + str(i) + ".list"
            s = Template( "-overwrite -in:file:native $native -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
            arguments = s.substitute( modellist=modelslist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)
            if ( restore == True ): 
                arguments = arguments + " -restore_talaris_behavior"

            # Generate a condor submission file and submit the job to Jazz
            condor_case_name = case + "_models_" + str(i)
            print "Submitting decoy-discrimination test case from Yarov-Yaravoy set:", condor_case_name
            if ( cluster_type == "SLURM" ): 
                write_and_submit_slurm_batch_script( outdir, condor_case_name, executable, arguments )
            else: 
                write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )


    ### Test Set #2: Dutagaci High Resolution Decoys (Molecular dynamics generated)
    list_of_test02_cases = inputs + "/dutagaci-set/decoy_sets.list"
    with open( list_of_test02_cases, 'rb' ) as f:
        test02_cases = f.readlines()
    test02_cases = [ x.strip() for x in test02_cases ]

    outdir_test02 = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/dutagaci-set"
    os.system( "mkdir " + outdir_test02 )
    os.chdir( outdir_test02 )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test02_cases:

        outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/dutagaci-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        native = inputs + "/dutagaci-set/" + case + "/" + case + "_native.pdb"
        spanfile = inputs + "/dutagaci-set/" + case + "/" + case + ".span"
        modelslist = inputs + "/dutagaci-set/" + case + "/decoys.list"
        s = Template( " -overwrite -in:file:native $native -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
        arguments = s.substitute( modellist=modelslist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"

        # Generate a condor submission file and submit the job to Jazz
        condor_case_name = case + "_dutagaci"
        print "Submitting decoy-discrimination test case from Dutagaci set:", condor_case_name
        if ( cluster_type == "SLURM" ): 
            write_and_submit_slurm_batch_script( outdir, condor_case_name, executable, arguments )
        else: 
            write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_t01 --which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option('--stable', '-s', 
        action="store", 
        help="For reference runs, use the stable Rosetta master branch", )

    parser.add_option('--cluster_type', '-t', 
        action="store", 
        help="Specify option to run jobs on a slurm or condor system", )

    parser.add_option( '--which_tests', '-w', 
        action="store", 
        help="Specify which test groups to run. Options are: ddG, landscape, prediction", )

    parser.add_option( '--restore_talaris', '-r', 
        action="store", 
        help="Restore talaris behavior using tthe flag -restore_talaris_behavior for reference runs"
        )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options are provided and contents are valid
    if ( not Options.energy_fxn or not Options.stable or not Options.cluster_type ): 
        print "Missing required options --energy_fxn, --stable, and/or --cluster_type! Exiting..."
        sys.exit()

    if ( ( not Options.cluster_type == "SLURM" ) and ( not Options.cluster_type == "CONDOR" ) ): 
        print "Invalid option for --cluster_type. Currently only supporting SLURM or CONDOR" 
        sys.exit()

    test_types = []
    if ( not Options.which_tests ): 
        test_types.append( "ddG" )
        test_types.append( "landscape" )
        test_types.append( "prediction" )
    else: 
        which_types = Options.which_tests.split(",")
        for t in which_types: 
            if ( t == "ddG" or t == "landscape" or t == "prediction" ): 
                test_types.append( t )
            else: 
                print "Invalid test type", t, "Exiting..."
                sys.exit()

    # Choose rosetta path based on input options
    rosetta_exe_path = ""
    if ( Options.stable == "true" ): 
        rosetta_exe_path = rosettadir_stable
    else: 
        rosetta_exe_path = rosettadir

    # Set option for restoring talaris behaviorrun_fixed
    restore = False
    if ( Options.restore_talaris == "true" ): 
        restore = True

    # If it doesn't exist, make the output data directory
    datadir = benchmark + "data/"
    if ( not os.path.isdir(datadir) ): 
        os.system( "mkdir " + datadir )
    os.system( "cd " + datadir )

    # Make an output data directory for the specific energy function
    outdir = benchmark + "data/" + Options.energy_fxn
    if ( not os.path.isdir(outdir) ): 
        os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    # Run energy landscape testing group
    if ( "landscape" in test_types ): 
    
        # Energy landscape test for single TM peptides found in nature
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.1_monomer_landscape", "helices.list", "xml/test_1.1_monomer_landscape.xml", restore )

        # Energy landscape test for aromatic-capped peptides
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.2_aro_landscape", "aro_helices.list", "xml/test_1.2_aro_landscape.xml", restore, "true" )

        # Energy landscape test for leucine-lysine peptides
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.3_lk_landscape", "lk_peptides.list", "xml/test_1.3_lk_landscape.xml", restore, "true" )

    # Run ddG calculations
    if ( "ddG" in test_types ): 

        # ddG of insertion landscape calculation for Ulmschneider set
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test-ddG-of-insertion", "insertion_peptide.dat", "xml/test_2.2_ddG_insertion_landscape.xml", restore, "true" )

        # ddG of insertion landscape calculation for pH dependent set - generate at pH = 4
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test-pH-dependent-insertion", "pH-inserted-helices.list", "xml/test_2.3_pH_landscape.xml", restore, "true", "4" )
    
        # ddG of insertion landscape calculation for pH dependent set - generate at pH = 7
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test-pH-dependent-insertion", "pH-inserted-helices.list", "xml/test_2.3_pH_landscape.xml", restore, "true", "7" )

        # ddG of mutation calculation for Moon & Fleming Set
        run_ddG_of_mutation_calc( Options.energy_fxn, "OmpLA/OmpLA_Moon_Fleming_set.dat", "OmpLA_Moon_Fleming_set", restore )

        # ddG of mutation calculation for McDonald & Fleming Set
        run_ddG_of_mutation_calc( Options.energy_fxn, "OmpLA_aro/OmpLA_aro_McDonald_Fleming_set.dat", "OmpLA_aro_McDonald_Fleming_set", restore )

        # ddG of mutation calculation for Marx & Fleming set
        run_ddG_of_mutation_calc( Options.energy_fxn, "PagP/PagP_Marx_Fleming_set.dat", "PagP_Marx_Fleming_set", restore )

    # Run prediction calculations
    if ( "prediction" in test_types ): 

        # Fixed backbone design calculation for sequence recovery test
        run_fixed_backbone_design_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, restore )

        # Docking calculation for small homodimer set (Lomize et al. 2017)
        run_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "small-homodimer-set", restore )

        # Docking calculation for large homodimer set (Alford & Koehler Leman 2015)
        run_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-homodimer-set", restore )

        # This doesn't have a label on it - so I'm wondering if this is where I had left off... 
        run_decoy_discrimination_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, restore )

if __name__ == "__main__" : main(sys.argv)
