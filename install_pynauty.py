import platform
import pip
import os
import sys
import tarfile
import shutil
import warnings
import argparse
import traceback

parser = argparse.ArgumentParser(description='A program to install pynauty')
parser.add_argument('--venv', help='define if inside a virtual environment', action="store_true")
parser.add_argument('--cmd_env', help='define if python find the correct environment for building [windows]', action="store_true")
args = parser.parse_args()
venv = bool(args.venv)
cmd_env = bool(args.cmd_env)

# Define operating system
is_linux, is_windows = False, False
if platform.system() == 'Linux':
    is_linux = True
elif platform.system() == 'Windows':
    # Warning: For a windows operating system MinGW must be installed.
    is_windows = True
else:
    sys.stderr.write('Unsupported os for this library')
    sys.exit(1)

# Running python version
os.system('echo python ' + str(platform.python_version()))
os.system('echo pip ' + str(pip.__version__))
python_executable_address = str(sys.executable)

# Download library
# for windows curl must be installed
if is_windows:
    os.system('curl https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz -o pynauty-0.6.0.tar.gz')
elif is_linux:
    os.system('wget https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz')
    
# Decompress
os.system('tar -xf pynauty-0.6.0.tar.gz')
os.remove('pynauty-0.6.0.tar.gz')

# Move inside pynauty
os.chdir("pynauty-0.6.0")

# Download C library
if is_windows:
    os.system('curl http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz -o nauty26r10.tar.gz')
elif is_linux:
    os.system('wget http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz')

# Decompress library
os.system('echo "decompressing: nauty26r10..."')
os.system('tar -xf nauty26r10.tar.gz')

if (is_windows):
    # rename folder
    os.system('echo "rename nauty"')
    os.system('rename nauty26r10 nauty')
    regular_exp = '\"\\\\?\\'+os.path.abspath('nauty')+'\\This_is_nauty26r10.\"'
    
    # Delete problematic file ending with '.'
    os.system('echo "delete config"')
    os.system('del '+regular_exp)

    # build pynauty
    os.system('make nauty-objects')
    if cmd_env:
        # Set environment variables
        pv = platform.python_version()
        old_pv = os.environ.get('PYTHON_VERSION','')
        os.environ['PYTHON_VERSION'] = pv

        # execute in env
        os.system('cmd /E:ON /V:ON /C .\\ci_scripts\\appveyor\\run_with_env.cmd ' + python_executable_address + ' setup.py build')
        
        # Restore
        if old_pv != '':
            os.environ['PYTHON_VERSION'] = old_pv
    else:
        os.system(python_executable_address + ' setup.py build')

if (is_linux):
    # rename folder
    os.system('mv nauty26r10 nauty')
    # build pynauty
    os.system('make nauty-objects')
    os.system(python_executable_address + ' setup.py build')

# Install pynauty
if venv:
    pip.main(['install','--upgrade','.'])
else:
    pip.main(['install','--user','--upgrade','.'])

# exit directory and delete
os.chdir("..")
shutil.rmtree('pynauty-0.6.0')
