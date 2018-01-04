import platform
import pip
import os
import sys
import tarfile
import shutil
import warnings
import argparse
import traceback

# Create an argument parser for the installer of pynauty
parser = argparse.ArgumentParser(description='A program to install pynauty')
parser.add_argument('--venv', help='define if inside a virtual environment', action="store_true")
meg = parser.add_mutually_exclusive_group()
meg.add_argument('--cmd_env', help='define if python needs to find the correct environment for building in visual studio [windows]', action="store_true")
meg.add_argument('--use_mingw', help='define if python needs to use the mingw compiler [windows]', action="store_true")
args = parser.parse_args()

# Define parameters useful for installation of pynauty as a consequence of argparse
venv = bool(args.venv)
cmd_env = bool(args.cmd_env)
use_mingw = bool(args.use_mingw)
compiler = ''
if use_mingw:
    compiler = ' --compiler=mingw32'

# Define operating system
is_linux, is_windows = False, False
if platform.system() == 'Linux':
    os.system('echo "Installing for Linux.."')
    is_linux = True
    if cmd_env:
        os.system('echo "Parameter --cmd_env is ignored.. Installation is happening inside a linus environment"')
    elif use_mingw:
        os.system('echo "Parameter --use_mingw is ignored.. Installation is happening inside a linus environment"')

elif platform.system() == 'Windows':
    # Warning: For a windows operating system MinGW must be installed.
    # Right now there is no support of pynauty for visual studio
    os.system('echo Installing for windows..')
    if cmd_env or not(use_mingw):
        os.system('echo Pynauty library is not supported for visual c++ compiler!')
        os.system('echo Switching to mingw32..')
    
    # Bypass parameter settings
    use_mingw = True
    cmd_env = False
    
    # Cygwin does not compile either.. Raises error:
    # C:\cygwin64\bin\gcc.exe -mcygwin -mdll -O -Wall -Inauty -Isrc -IC:\Python34\include -IC:\Python34\include -c src/nautywrap.c -o build\temp.win-amd64-3.4\Release\src\nautywrap.o -O4 -fPIC
    # gcc: error: unrecognized command line option '-mcygwin'
    compiler = ' --compiler=mingw32'

    if use_mingw:
        os.system('echo Importing a mingw32 build library on a python34 system')
        os.system('echo the only valid system for intallation of grakel on Windows (because of cvxopt dependencies..)')
        os.system('echo allocates a 4 GB memory space on import..')
        os.system('echo Even if that is acceptable system-wise, the execution speed will propably decrease..')
        os.system('echo This is not a good build..')
        os.system('echo For this issue see here:')
        os.system('echo https://github.com/ContinuumIO/anaconda-issues/issues/271#issue-58658137')
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
    os.system('echo Making nauty object files ..')
    os.system('make nauty-objects')
    os.system('echo Object Files made ..')

    if cmd_env:
        # Set environment variables
        pv = platform.python_version()
        old_pv = os.environ.get('PYTHON_VERSION','')
        os.environ['PYTHON_VERSION'] = pv

        os.system('echo Building inside environment')
        
        # Execute with the correct environment
        # problem!! nauty.h is designed for gcc compilation
        os.system('cmd /E:ON /V:ON /C ..\\ci_scripts\\appveyor\\run_with_env.cmd ' + python_executable_address + ' setup.py build')
        
        # Restore
        if old_pv != '':
            os.environ['PYTHON_VERSION'] = old_pv
    else:
        os.system(python_executable_address + ' setup.py build' + compiler)

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
