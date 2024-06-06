#!/usr/bin/env python3
import os
import glob
import sys
import numpy as np
import argparse
import logging
import time
import subprocess
import json
import copy


pt = logging.warning
verb = logging.info

# CL Argument Parser
parser = argparse.ArgumentParser()
parser.add_argument('--verbose','-v', action='store_true', help='Verbose output, default off')
parser.add_argument('--subdirs','-s', action='store_false', help='Run in all direct subdirectories, default on')
parser.add_argument('--overwrite','-o', action='store_true', help='Overwrite existing slurm scripts, default off')
parser.add_argument('--input','-i', action='store_true', help='SLURM header manual inputs. Default=Off')
parser.add_argument('--inl','--pbs','-p', action='store_true', help='PBS script for INL. Default=Off')
parser.add_argument('--run','-r', action='store_true', help='SBATCH Submit all the slurm scripts also. Default=Off')
parser.add_argument('--json','-j', action='store_true', help='Read SLURM header values from SLURM.json, Default=Off')
parser.add_argument('--dest','-d', action='store_true', help='''DISABLED: Change the file_base in the MOOSE inputs to
                    match the /blue equivalent with the input.i name as the output name still, Default=Off''')

# Slurm Header Args
parser.add_argument('--dir-names', action='store_false', help='''SLURM Job Name set to directory names, defaults to true.
                    If false will use .i file names''')
parser.add_argument('--nodes','-n', type=int, help='SLURM number of nodes to use. Default=1')
parser.add_argument('--tasks', type=int, help='SLURM number of (mpi?) tasks. Default=30')
parser.add_argument('--cpus-per-task', type=int, help='SLURM number of cpus per task to use. Default=1')
parser.add_argument('--mem-per-cpu', type=str, help='''SLURM memory per cpu. The default is fairy safe.
                    Default=7GB''')
parser.add_argument('--partition', type=str, help='SLURM partition to use. Default=NONE')
parser.add_argument('--burst','-b', action='store_true', help='SLURM run as a burst job, default off')
parser.add_argument('--time','-t', type=int, help='''SLURM number of hours to run, on the default partitions
                    of [hpg-default, hpg2-compute, bigmem] burst limit is 96, while regular limit is 744 (31 days).
                    Default=72 hours''')
parser.add_argument('--args','-c', type=str, help='Extra CL arguments at the end. Default=NONE')
cl_args = parser.parse_args()

# Defaults for the variables
class default_vals:
    dir_names = True
    nodes = 1
    tasks = 30
    cpus_per_task = 1
    mem_per_cpu = '7GB'
    time = 72
    partition = None
    burst = False
# save the current input values of cl_args for reference
input_vals = copy.copy(cl_args)


# Toggle verbose
if cl_args.verbose == True:
    logging.basicConfig(level=logging.INFO,format='%(message)s')
elif cl_args.verbose == False:
    logging.basicConfig(level=logging.WARNING,format='%(message)s')
verb('Verbose Logging Enabled')
verb(cl_args)


# Some Warnings
if cl_args.subdirs:
    verb('By default this script runs in 1 level of subdirectories')
else:
    verb('Only running in current directory, not in subdirectories.')
if cl_args.overwrite:
    pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Overwriting existing .sh files')
if cl_args.run:
    pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Running all relevant .sh files')
if cl_args.json:
    pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' command line flagged values will'
       +' be overwritten by those in the json file')
if cl_args.dest:
    pt('\x1b[31;1m'+'ERROR:'+'\x1b[0m'+' NOT SET UP: Changing file_base in the input files to'
       +' direct output to another loaction (/blue on HPG likely)')
pt(' ')


# Pathing
# Find absolute path to phase_field-opt executable
# pf_opt = os.path.expanduser('~/projects/moose/modules/phase_field/phase_field-opt')
script_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
pf_opt = script_directory.rsplit('/',3)[0]+'/phase_field-opt'
if os.path.exists(pf_opt):
    verb('Found Phase_Field-opt')
else:
    pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Cannot find phase_field-opt, may cause issues writing SLURM script')
    confirm = input('Would you like to continue anyway? y/[n]: ')
    if 'y' not in confirm.lower():
        pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Exiting Code')
        sys.exit()
if cl_args.json:
    json_path = script_directory+'/slurm_header_template.json'
    if os.path.exists(json_path):
        verb('Found slurm_header_template.json')
        write_json = False
    else:
        pt('\x1b[31;1m'+'ERROR:'+'\x1b[0m'+' Cannot find slurm_header_template.json,'+
           ' will cause issues writing SLURM script')
        pt('Script path should be: '+json_path)
        confirm = input('Would you like to write the current output to a new json file? [y]/n: ')
        if 'n' not in confirm.lower():
            pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Writing a new slurm_header_template.json')
            write_json = True
        else:
            confirm = input('Would you like to continue anyway? y/[n]: ')
            if 'y' not in confirm.lower():
                pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Exiting Code')
                sys.exit()
            else:
                pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Bad idea, good luck...')
                write_json = False
verb(' ')


# ███████╗██╗   ██╗███╗   ██╗ ██████╗████████╗██╗ ██████╗ ███╗   ██╗███████╗
# ██╔════╝██║   ██║████╗  ██║██╔════╝╚══██╔══╝██║██╔═══██╗████╗  ██║██╔════╝
# █████╗  ██║   ██║██╔██╗ ██║██║        ██║   ██║██║   ██║██╔██╗ ██║███████╗
# ██╔══╝  ██║   ██║██║╚██╗██║██║        ██║   ██║██║   ██║██║╚██╗██║╚════██║
# ██║     ╚██████╔╝██║ ╚████║╚██████╗   ██║   ██║╚██████╔╝██║ ╚████║███████║
# ╚═╝      ╚═════╝ ╚═╝  ╚═══╝ ╚═════╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝

# Check current directory for .sh script, true if it exists
def checkForSlurm():
    if glob.glob("*.sh"):
        if cl_args.overwrite:
            verb('  Overwriting .sh script')
            return False
        else:
            verb('  Found .sh script')
            return True
    else:
        return False


# Check current directory for .i and if it exists return the name - .i
def checkForInput():
    if glob.glob("*.i"):
        for file in glob.glob("*.i"):
            verb('  Input File is: ' + file)
            inputName = file.rsplit('.',1)[0]
        return True, inputName
    else:
        verb('  No Input File')
        return False, False


# Optional manual input mode for the first one
def manualInput():
    pt('Inputting SLURM header values manually')
    pt('  An empty entry will use the default or current value')
    # Job Name
    tempVal = 'Input Name'
    if cl_args.dir_names:
        tempVal = 'Directory Name'
    slname = input('Job Name (dir or input) (Currently '+tempVal+'): ')
    if 'd' in slname:
        verb('  Using Directory Names')
        cl_args.dir_names = True
    elif slname == '':
        verb('  Using current value of '+tempVal)
    else:
        verb('  Using Input.i Names')
        cl_args.dir_names = False
    # Nodes
    nodes = input('Number of nodes (currently '+str(cl_args.nodes)+'): ')
    if nodes == '':
        verb('  Using current value of '+str(cl_args.nodes)+' nodes')
    else:
        cl_args.nodes = nodes
    # Tasks/CPUs per node
    ncpu = input('Number of CPUs per node (currently '+str(cl_args.tasks)+'): ')
    if ncpu == '':
        verb('  Using current value of '+str(cl_args.tasks)+' CPUs per node')
    else:
        cl_args.tasks = ncpu
        cl_args.cpus_per_task = 1
    # Memory
    mem = input('Memory per CPU (currently '+str(cl_args.mem_per_cpu)+'): ')
    if mem == '':
        verb('  Using current value of '+str(cl_args.mem_per_cpu)+' per CPU')
    else:
        cl_args.mem_per_cpu = mem
    # Partition
    partition = input('Partition (optional, None for unspecified, currently '+str(cl_args.partition)+'): ')
    if partition == '':
        verb('  Using current value of '+str(cl_args.partition))
    else:
        cl_args.partition = partition
    # Time
    time = input('Time to run in hours (or days in format d-hh, currently '+str(cl_args.time)+'): ')
    if time == '':
        verb('  Using current value of '+str(cl_args.time)+':00:00')
    else:
        cl_args.time = time
    # Burst
    burst = input('Burst (TF, currently '+str(cl_args.burst)+'): ')
    if burst == '':
        verb('  Using current value of '+str(cl_args.burst))
    else:
        cl_args.burst = burst
    slurmHeaderPreview(True)
    pt(' ')
    return


# Write the actual slurm script
def slurmWrite(cwd,inputName):
    slurmList = []
    if cl_args.inl == False:
        # Slurm style
        # Building the header
        slurmList.append('#!/bin/bash')
        slurmList.append('')
        # Job Name: Directory or Input File
        if cl_args.dir_names:
            slurmList.append('#SBATCH --job-name='+cwd.rsplit('/',1)[1])
        elif cl_args.dir_names == False:
            slurmList.append('#SBATCH --job-name='+inputName)
        else:
            raise ValueError("--dir-names not specified True or False")
        # Nodes
        slurmList.append('#SBATCH --nodes='+str(cl_args.nodes))
        # Tasks per node
        slurmList.append('#SBATCH --ntasks-per-node='+str(cl_args.tasks))
        # CPUs per task
        slurmList.append('#SBATCH --cpus-per-task='+str(cl_args.cpus_per_task))
        # Memory per CPU
        slurmList.append('#SBATCH --mem-per-cpu='+cl_args.mem_per_cpu)
        # Distribution line
        slurmList.append('#SBATCH --distribution=cyclic:cyclic')
        # A nice empty line for cleanliness
        slurmList.append('')
        # Partition if specified
        if cl_args.partition==None or cl_args.partition=="None":
            verb('    No partition specified')
        else:
            slurmList.append('#SBATCH --partition='+cl_args.partition)
        # Time to run in hours
        slurmList.append('#SBATCH --time='+str(cl_args.time)+':00:00')
        # Terminal save name
        slurmList.append('#SBATCH --output=moose_console_%j.out')
        # Email
        slurmList.append('#SBATCH --mail-user=bbattas@ufl.edu')
        slurmList.append('#SBATCH --mail-type=BEGIN,END,FAIL')
        # Account to run on and burst or not
        slurmList.append('#SBATCH --account=michael.tonks')
        if cl_args.burst:
            verb('    Specifying burst allocation')
            slurmList.append('#SBATCH --qos=michael.tonks-b')
        # Current Exclude List (2/1/24)
        # slurmList.append('#SBATCH --exclude=c0702a-s28,c0702a-s29,c0703a-s18,c0706a-s7,c0709a-s21,c0710a-s28,c0713a-s18,c0713a-s19')

        # On to the actual job to submit
        # Define Locations
        slurmList.append('')
        slurmList.append('echo ${SLURM_JOB_NODELIST}')
        slurmList.append('')
        slurmList.append('MOOSE='+pf_opt)
        slurmList.append('OUTPUT='+cwd)
        # Module loading
        slurmList.append('')
        slurmList.append('export CC=mpicc CXX=mpicxx FC=mpif90 F90=mpif90 F77=mpif77')
        slurmList.append('module purge')
        slurmList.append('module load ufrc mkl/2023.2.0 gcc/12.2.0 openmpi/4.1.5 python/3.11 cmake/3.26.4')
        # slurmList.append('module load conda')
        # Mamba Build
        # slurmList.append('source ~/.bashrc')
        # slurmList.append('mamba activate moose')
        # Actually go to the output and run the shit
        slurmList.append('')
        slurmList.append('cd $OUTPUT')
        if cl_args.args == None:
            slurmList.append('srun --mpi=pmix_v3 $MOOSE -i $OUTPUT/'+inputName+'.i')
        else:
            slurmList.append('srun --mpi=pmix_v3 $MOOSE -i $OUTPUT/'+inputName+'.i '+str(cl_args.args))

        # Output the slurm script
        # verb(slurmList)
        slurmName = 'slurm_' + inputName + '.sh'
        verb('    Saving slurm script: '+slurmName)
        with open(slurmName, mode='w') as slurmFile:
            slurmFile.write('\n'.join(slurmList))
    else:
        # PBS style
        # Building the header
        slurmList.append('#!/bin/bash')
        # Job Name: Directory or Input File
        if cl_args.dir_names:
            slurmList.append('#PBS -N '+cwd.rsplit('/',1)[1])
        elif cl_args.dir_names == False:
            slurmList.append('#PBS -N '+inputName)
        else:
            raise ValueError("--dir-names not specified True or False")
        # Nodes
        slurmList.append('#PBS -l select='+str(cl_args.nodes)+':ncpus=48:mpiprocs=48')
        # Time to run in hours
        slurmList.append('#PBS -l walltime='+str(cl_args.time)+':0:0')
        slurmList.append('#PBS -k doe')
        slurmList.append('#PBS -o terminal_out')
        slurmList.append('#PBS -j oe')
        slurmList.append('#PBS -P neams')
        slurmList.append('#PBS -m abe')
        slurmList.append('#PBS -M bbattas@ufl.edu')
        slurmList.append(' ')
        slurmList.append('cat $PBS_NODEFILE')
        slurmList.append('source /etc/profile.d/modules.sh')
        slurmList.append('module load use.moose moose-dev')
        slurmList.append(' ')
        slurmList.append('cd $PBS_O_WORKDIR')
        slurmList.append(' ')
        # slurmList.append('mpirun /home/$USER/projects/moose/modules/phase_field/phase_field-opt -i '+inputName+'.i')
        slurmList.append('mpirun /scratch/$USER/moose/modules/phase_field/phase_field-opt -i '+inputName+'.i')
        # Output the PBS script
        # verb(slurmList)
        slurmName = 'pbs_' + inputName + '.sh'
        verb('    Saving slurm script: '+slurmName)
        with open(slurmName, mode='w') as slurmFile:
            slurmFile.write('\n'.join(slurmList))


# Preview of the header content of the SLURM script printed to terminal
# copied from the slurmWrite(), so if the header there changes it needs to here too
def slurmHeaderPreview(interactive):
    if cl_args.inl == True:
        pt('No Preview for PBS currently')
        return
    if interactive:
        pt(' ')
        pt('\x1B[1m'+'Header Preview: \x1b[0m')
        pt(' ')
    slurmList = []
    # Building the header
    slurmList.append('#!/bin/bash')
    slurmList.append('')
    # Job Name: Directory or Input File
    if cl_args.dir_names:
        slurmList.append('#SBATCH --job-name='+'[CWD]')
    elif cl_args.dir_names == False:
        slurmList.append('#SBATCH --job-name='+'[inputName]')
    else:
        raise ValueError("--dir-names not specified True or False")
    # Nodes
    slurmList.append('#SBATCH --nodes='+str(cl_args.nodes))
    # Tasks per node
    slurmList.append('#SBATCH --ntasks-per-node='+str(cl_args.tasks))
    # CPUs per task
    slurmList.append('#SBATCH --cpus-per-task='+str(cl_args.cpus_per_task))
    # Memory per CPU
    slurmList.append('#SBATCH --mem-per-cpu='+cl_args.mem_per_cpu)
    # Distribution line
    slurmList.append('#SBATCH --distribution=cyclic:cyclic')
    # A nice empty line for cleanliness
    slurmList.append('')
    # Partition if specified
    if cl_args.partition==None or cl_args.partition=="None":
        logging.debug('    No partition specified')
    else:
        slurmList.append('#SBATCH --partition='+cl_args.partition)
    # Time to run in hours
    slurmList.append('#SBATCH --time='+str(cl_args.time)+':00:00')
    # Terminal save name
    slurmList.append('#SBATCH --output=moose_console_%j.out')
    # Email
    slurmList.append('#SBATCH --mail-user=bbattas@ufl.edu')
    slurmList.append('#SBATCH --mail-type=BEGIN,END,FAIL')
    # Account to run on and burst or not
    slurmList.append('#SBATCH --account=michael.tonks')
    if cl_args.burst:
        # verb('    Specifying burst allocation')
        slurmList.append('#SBATCH --qos=michael.tonks-b')
    # # Current Exclude List (2/1/24)
    # slurmList.append('#SBATCH --exclude=c0702a-s28,c0702a-s29,c0703a-s18,c0706a-s7,c0709a-s21,c0710a-s28,c0713a-s18,c0713a-s19')
    pt('\n'.join(slurmList))
    if interactive:
        pt (' ')
        confirm = input('Continue using this as the header? [y]/n: ')
        if 'n' in confirm.lower():
            pt('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+' Exiting Code')
            sys.exit()


# Run the slurm script in the current directory
def runSlurm():
    if cl_args.inl == False:
        verb('  Submitting slurm script')
        # Find slurm script name
        for file in glob.glob("*.sh"):
            slurmScriptName = file
        command = ['sbatch',slurmScriptName]
        verb('    Command being run is: ' + str(command))
        # Submit the slurm script
        subprocess.run(command)
    else:
        verb('  Submitting PBS script')
        # Find slurm script name
        for file in glob.glob("*.sh"):
            slurmScriptName = file
        command = ['qsub',slurmScriptName]
        verb('    Command being run is: ' + str(command))
        # Submit the slurm script
        subprocess.run(command)


def setCL_valuesOverride(defaults_tf):
    verb('Setting cl_args based on CL overrides')
    # Set defaults first (if not loading json)
    if defaults_tf is True:
        cl_args.dir_names = default_vals.dir_names
        cl_args.nodes = default_vals.nodes
        cl_args.tasks = default_vals.tasks
        cl_args.cpus_per_task = default_vals.cpus_per_task
        cl_args.mem_per_cpu = default_vals.mem_per_cpu
        cl_args.time = default_vals.time
        cl_args.partition = default_vals.partition
        cl_args.burst = default_vals.burst
    # If the cl_args input values not None (or default for T/F)
    if input_vals.dir_names is False:
        cl_args.dir_names = input_vals.dir_names
    if input_vals.nodes is not None:
        cl_args.nodes = input_vals.nodes
    if input_vals.tasks is not None:
        cl_args.tasks = input_vals.tasks
    if input_vals.cpus_per_task is not None:
        cl_args.cpus_per_task = input_vals.cpus_per_task
    if input_vals.mem_per_cpu is not None:
        cl_args.mem_per_cpu = input_vals.mem_per_cpu
    if input_vals.time is not None:
        cl_args.time = input_vals.time
    if input_vals.partition is not None:
        cl_args.partition = input_vals.partition
    if input_vals.burst is True:
        cl_args.burst = input_vals.burst



# Read a json file 'slurm_header_template.json' in scripts
# folder under pf/examples/sintering/ to define the different
# slurm script header parameters
def jsonToSlurmVariables():
    # If --json and file exists
    if cl_args.json and write_json==False:
        verb('Reading from json')
        # verb('  '+json_path)
        # Read the json file
        with open(json_path) as json_file:
            dict = json.load(json_file)
        # Redefine cl_args based on the json values
        cl_args.dir_names = dict['job_name_using_dir_names']
        cl_args.nodes = dict['nodes']
        cl_args.tasks = dict['tasks_per_node']
        cl_args.cpus_per_task = dict['cpus_per_task']
        cl_args.mem_per_cpu = dict['mem_per_cpu']
        cl_args.time = dict['hours']
        cl_args.partition = dict['partition']
        cl_args.burst = dict['burst']
        setCL_valuesOverride(False)
        slurmHeaderPreview(True)
        # verb(' ')
        # verb(cl_args)
        # verb(' ')
    # If --json and file does not exist
    elif cl_args.json and write_json:
        verb('Writing a new json file at: ')
        verb('  '+json_path)
        if cl_args.input:
            manualInput()
        setCL_valuesOverride(True)
        dict = {}
        dict['job_name_using_dir_names'] = cl_args.dir_names
        dict['nodes'] = cl_args.nodes
        dict['tasks_per_node'] = cl_args.tasks
        dict['cpus_per_task'] = cl_args.cpus_per_task
        dict['mem_per_cpu'] = cl_args.mem_per_cpu
        dict['hours'] = cl_args.time
        dict['partition'] = cl_args.partition
        dict['burst'] = cl_args.burst

        with open(json_path, 'w') as fp:
            json.dump(dict, fp)
        slurmHeaderPreview(True)

    else:
        logging.debug('Skipping json as the -json flag = '+str(cl_args.json))
        verb('No json used, setting cl_args from defaults')
        setCL_valuesOverride(True)
        slurmHeaderPreview(True)



# ███████╗██╗  ██╗███████╗ ██████╗██╗   ██╗████████╗███████╗
# ██╔════╝╚██╗██╔╝██╔════╝██╔════╝██║   ██║╚══██╔══╝██╔════╝
# █████╗   ╚███╔╝ █████╗  ██║     ██║   ██║   ██║   █████╗
# ██╔══╝   ██╔██╗ ██╔══╝  ██║     ██║   ██║   ██║   ██╔══╝
# ███████╗██╔╝ ██╗███████╗╚██████╗╚██████╔╝   ██║   ███████╗
# ╚══════╝╚═╝  ╚═╝╚══════╝ ╚═════╝ ╚═════╝    ╚═╝   ╚══════╝


# Decision tree for manual input and json to overwrite cl_args
if cl_args.input and cl_args.json:
    if write_json:
        # Writing a new json with manual input
        jsonToSlurmVariables()
    else:
        # Read the json file at pf/examples/sintering/scripts
        jsonToSlurmVariables()
        pt(' ')
        pt('Manual overwriting the slurm values:  ')
        # Now manual overwrite of json input
        manualInput()
elif cl_args.input and not cl_args.json:
    manualInput()
elif cl_args.json and not cl_args.input:
    jsonToSlurmVariables()
elif not cl_args.input and not cl_args.json:
    verb('No json or manual input for SLURM header- using defaults with command line overrides')
    jsonToSlurmVariables()
    verb(' ')
else:
    raise ValueError('\x1b[31;1m'+'ERROR:'+'\x1b[0m'+' json and manual input logic tree failure')



# All subdirectories 1 level deep
if cl_args.subdirs:
    verb('Cycling through all subdirectories.')
    # Cycle through the directories
    if cl_args.subdirs:
        working_dir = os.getcwd()
        verb('Directories:')
        verb(glob.glob("*/"))
        for dir in os.listdir():
            if not dir.startswith('.'):
                cwd = working_dir + "/" + dir
                os.chdir(cwd)
                verb('In directory: ' + os.getcwd())
                # Check if there is a .i file
                input = checkForInput()
                if input[0]:
                    # Check if we need to make a slurm script for the .i file
                    if checkForSlurm() == False:
                        # Make a new slurm script
                        slurmWrite(cwd,input[1])
                    if cl_args.run:
                        runSlurm()
                    verb(' ')
# No sibdirectories, only the current directory
else:
    verb('Running only in current directory')
    cwd = os.getcwd()
    verb('In directory: ' + cwd)
    # Check if there is a .i file
    input = checkForInput()
    if input[0]:
        # Check if we need to make a slurm script for the .i file
        if checkForSlurm() == False:
            # Make a new slurm script
            slurmWrite(cwd,input[1])
        if cl_args.run:
            runSlurm()
        verb(' ')

verb(' ')
verb('██████╗  ██████╗ ███╗   ██╗███████╗')
verb('██╔══██╗██╔═══██╗████╗  ██║██╔════╝')
verb('██║  ██║██║   ██║██╔██╗ ██║█████╗  ')
verb('██║  ██║██║   ██║██║╚██╗██║██╔══╝  ')
verb('██████╔╝╚██████╔╝██║ ╚████║███████╗')
verb('╚═════╝  ╚═════╝ ╚═╝  ╚═══╝╚══════╝')
verb('''
⢀⣠⣾⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⠀⣠⣤⣶⣶
⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠀⠀⠀⢰⣿⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣧⣀⣀⣾⣿⣿⣿⣿
⣿⣿⣿⣿⣿⡏⠉⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⡿⣿
⣿⣿⣿⣿⣿⣿⠀⠀⠀⠈⠛⢿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⠿⠛⠉⠁⠀⣿
⣿⣿⣿⣿⣿⣿⣧⡀⠀⠀⠀⠀⠙⠿⠿⠿⠻⠿⠿⠟⠿⠛⠉⠀⠀⠀⠀⠀⣸⣿
⣿⣿⣿⣿⣿⣿⣿⣷⣄⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⣴⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⣿⠏⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠠⣴⣿⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⡟⠀⠀⢰⣹⡆⠀⠀⠀⠀⠀⠀⣭⣷⠀⠀⠀⠸⣿⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀⠈⠉⠀⠀⠤⠄⠀⠀⠀⠉⠁⠀⠀⠀⠀⢿⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⢾⣿⣷⠀⠀⠀⠀⡠⠤⢄⠀⠀⠀⠠⣿⣿⣷⠀⢸⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⡀⠉⠀⠀⠀⠀⠀⢄⠀⢀⠀⠀⠀⠀⠉⠉⠁⠀⠀⣿⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⣧⠀⠀⠀⠀⠀⠀⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢹⣿⣿
⣿⣿⣿⣿⣿⣿⣿⣿⣿⠃⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⣿⣿''')

verb(' ')


verb("Finished Exectuing Python Script")
quit()
quit()
