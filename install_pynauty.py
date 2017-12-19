import platform
import pip
import os
import sys
import tarfile
import shutil

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
pip.main(['install', '--upgrade', 'wget'])
if (is_windows):
    os.system('mingw64-get install make')

import wget
# Download library
wget.download('https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz')

# Decompress
tar = tarfile.open("pynauty-0.6.0.tar.gz", "r:gz")
tar.extractall()
tar.close()
os.remove('pynauty-0.6.0.tar.gz')

# Move inside pynauty
os.chdir("pynauty-0.6.0")

# Download C library
wget.download('http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz')

# Decompress library
tar = tarfile.open("nauty26r10.tar.gz", "r:gz")
tar.extractall()
tar.close()

# make a soft link
os.symlink('nauty26r10', 'nauty')


# build pynauty
os.system('make pynauty')

# define if insidd virtual-env and install
if hasattr(sys, 'real_prefix'):
    os.system('make virtenv-ins')
else:
    os.system('make tests')
    os.system('make user-ins')

# exit directory and delete
os.chdir("..")
shutil.rmtree('pynauty-0.6.0')
