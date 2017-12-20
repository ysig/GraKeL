import platform
import pip
import os
import sys
import tarfile
import shutil
import warnings
import argparse

parser = argparse.ArgumentParser(description='A program to install pynauty')
parser.add_argument('--venv', help='define if inside a virtual environment', action="store_true")
args = parser.parse_args()
venv = bool(args.venv)

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

# Dependencies
pip.main(['install', '--upgrade', 'setuptools'])

# Download library
os.system('wget https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz')

# Decompress
os.system('tar -xf pynauty-0.6.0.tar.gz')
os.remove('pynauty-0.6.0.tar.gz')

# Move inside pynauty
os.chdir("pynauty-0.6.0")

# Download C library
# wget must be installed
os.system('wget http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz')

# Decompress library
os.system('tar -xf nauty26r10.tar.gz')

# rename folder
if (is_windows):
    os.system('rename nauty26r10 nauty')
    print('del \\\\?'+os.path.join(os.path.abspath("nauty"),"This_is_nauty_26r10."))
    os.system('del \\\\?'+os.path.join(os.path.abspath("nauty"),"This_is_nauty_26r10."))
    os.system('dir .')
    os.system('dir nauty')
if (is_linux):
    os.system('mv nauty26r10 nauty')

# build pynauty
if is_windows:
    os.system('make nauty-objects')
    os.system('python setup.py build --compiler=mingw32')
if is_linux:
    os.system('make pynauty')
    
print("PyNauty Build Succesfully!")
# define if inside virtual-env and install
if venv:
    os.system('pip install --upgrade .')
else:
    os.system('pip install --user --upgrade .')
# exit directory and delete
os.chdir("..")
shutil.rmtree('pynauty-0.6.0')
