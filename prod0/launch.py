import os
import subprocess

ISOL='TOPO_PTDEP'
input_dir = 'LocalInput/runs'
output_dir = '/sps/atlas/q/qbuat/BUMPY_BUMP_REL17ISO/'

def cmd(
    input, 
    output,
    isol='TOPO',
    driver='submit.sh',
    name='toto'):
    
    cmd_line = 'qsub -N {0} {1} {2} {3} {4}'.format(
        name, driver, isol,
        os.path.join(input_dir, input), 
        os.path.join(output_dir, output))
    return cmd_line


for file in os.listdir(input_dir):
    if not 'run' in file:
        continue
    # print file
    run = file.split('_')[-1].split('.')[0]
    output = 'datafile_' + ISOL + '_' + run + '.root'

    command = cmd(file, output, isol=ISOL, name=run)
    print command
    subprocess.call(command, shell=True)
    # break
    # print output
    
