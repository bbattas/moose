import subprocess
import os
import glob
import argparse


# Command Line Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--cpus','-n',type=int,default=2,
                    help='Number of CPUs to run on (Default=2)')
parser.add_argument('--dim','-d',type=int,default=2, choices=[2,3],
                    help='Dimensions for the calculation (Default=2)')
parser.add_argument('--time','-t', action='store_true',
                            help='Run time_file_make script, default off')
parser.add_argument('--skip','-s', action='store_true',
                            help='''In the time_file_make use the skip flag '''
                            '''to skip last file/timestep(s), default off''')
parser.add_argument('--max_xy','-x',type=float,required=True,
                            help='''Max x&y for the 1/4 hull assuming center of symmetry is '''
                            '''[max_xy,max_xy,max_z], Required, Default=None''')
parser.add_argument('--max_z','-z',type=float,default=None,
                            help='''Max z for the 1/4 hull assuming center of symmetry is '''
                            '''[max_xy,max_xy,max_z], Required for 3D, Default=None''')
cl_args = parser.parse_args()



# print('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+
#       ' Currently the dimensions for the inputs for the quarter hull math need to be manually '+
#       'entered in the TIGER script!  Make sure you did this before running each size of input!!!')

if cl_args.dim == 3:
    if cl_args.max_z == None:
        raise ValueError('For 3D, --max_z/-z is required!')

working_dir = os.getcwd()
#os.chdir(working_dir)
print(working_dir)
print("Directories to run in:")
#print(os.listdir())
# print(glob.glob("*/"))
dir_names = sorted(glob.glob('*/'))
# print(dir_names)
for n in dir_names:
    print(n)
#subprocess.call(["ls","-d","--","*/"])#'ls -d -- */'
print(" ")

# File Paths
full_time = os.path.abspath(os.path.expanduser('~/projects/TIGER/parallel_time_file_make.py'))
full_2d = os.path.abspath(os.path.expanduser('~/projects/TIGER/TEMP_parallel_2d_area_from_timeFile.py'))
full_3d = os.path.abspath(os.path.expanduser('~/projects/TIGER/TEMP_parallel_3d_volume_from_timeFile.py'))

# Build command(s)
time_command = ['python',full_time,str(cl_args.cpus)]
if cl_args.skip:
    time_command.append('skip')
area_command = ['python',full_2d,'-n',str(cl_args.cpus),'-x',str(cl_args.max_xy)]
vol_command = ['python',full_3d,'-n',str(cl_args.cpus),'-x',str(cl_args.max_xy),'-z',str(cl_args.max_z)]

print('Commands being run: ')
if cl_args.time:
    print(time_command)
if cl_args.dim == 2:
    print(area_command)
if cl_args.dim == 3:
    print(vol_command)
print(' ')
print(' ')

#for dir in os.listdir():
for dir in dir_names:
    if not dir.startswith('.') and not dir.startswith('pic'):# and not dir.startswith('3D'):
        os.chdir(working_dir + "/" + dir)
        print("In directory: ",os.getcwd())
        if cl_args.time:
            subprocess.run(time_command)
        if cl_args.dim == 2:
            subprocess.run(area_command)
        if cl_args.dim == 3:
            subprocess.run(vol_command)
        # subprocess.run(["python","/home/bbattas/projects/TIGER/parallel_time_file_make.py","30"])#,"skip"
        # subprocess.run(["python","/home/bbattas/projects/TIGER/TEMP_parallel_3d_volume_from_timeFile.py","30"])
        # subprocess.run(["python","/home/bbattas/projects/TIGER/TEMP_parallel_2d_area_from_timeFile.py","60"])
        # subprocess.run(["python","/home/bbattas/projects/TIGER/parallel_2d_area_from_phi_timeFile.py","60"])
        # subprocess.run(["python","/home/bbattas/projects/TIGER/3d_plane_from_timeFile.py"])
        os.chdir(working_dir)
print(" ")
print(" ")
print("Done")
print(" ")
print('\x1b[31;1m'+'WARNING:'+'\x1b[0m'+
      ' Currently the dimensions for the inputs for the quarter hull math need to be manually '+
      'specified in the command line flags for this script!  Remember there are two sizes of inputs...')
print(" ")
