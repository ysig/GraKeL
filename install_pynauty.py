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
    overwrite = True
    exist_cfg_bak = False
    config_lines = '[build]\ncompiler=mingw32\n[build_ext]\ncompiler=mingw32'

    # For non virtual env propmpt
    try: input = raw_input
    except NameError: pass
    if not venv:
        print('To setup pynauty the mingw32 compiler must be used.')
        print('The following lines must added or substituted to "python_install_dir"\\Lib\\distutils\\distutils.cfg')
        valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
        sys.stdout.write("Do you want to overwrite (if exists)?" + " [y/N]")
        choice = input().lower()
        if choice not in valid: 
            overwrite = False
        else:
            overwrite = not valid[choice]
    
    # Add build-compiler variables in config.cfg
    if overwrite:
        f = False
        if os.environ.get('PYTHON', '')=='':
            print('You must set a system variable PYTHON pointing to your python installed directory')
            f = True
        elif not os.path.exists(os.environ['PYTHON'] + '\\Lib\\distutils'):
            f = False
            sys.stderr.write('The current system variable PYTHON pointing to your python installed directory is invalid\nPython Directory must contain \\Lib\\distutils ..')
        while f:
            print('SET PYTHON=')
            choice = input()
            f = os.path.exists(str(choice) + '\\Lib\\distutils')
            if not f:
                sys.stderr.write(str(choice) + ' is not a valid Python directory\nPython Directory must contain \\Lib\\distutils ..')
            else:
                os.system('set PYTHON='+str(choice))
                os.environ['PYTHON'] = str(choice)
        config_path = str(os.environ['PYTHON']) + '\\Lib\\distutils\\distutils.cfg'
        if os.path.exists(config_path):
            warning.warn('File distutils.cfg exists..')
            os.system('copy ' + config_path + ' ' + config_path + '.bak')
            warning.warn('Overwriting .. [Backup saved on "distutils.cfg.bak"]')
            os.system('del ' + config_path)
            exist_cfg_bak = True
        open(config_path, 'wb').write('[build]\ncompiler=mingw32\n[build_ext]\ncompiler=mingw32')
        
if (is_linux):
    os.system('mv nauty26r10 nauty')

# build pynauty
os.system('make pynauty')

# define if insidd virtual-env and install
if venv:
    os.system('make virtenv-ins')
else:
    os.system('make tests')
    os.system('make user-ins')

# exit directory and delete
os.chdir("..")
shutil.rmtree('pynauty-0.6.0')
if is_windows and exists_cfg_bak:
    warning.warn('Restoring: distutils.cfg.bak -> distutils.cfg')
    os.system('del ' + config_path)
    os.system('rename ' + config_path + '.bak ' + config_path)
