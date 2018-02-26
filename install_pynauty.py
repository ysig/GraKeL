"""A python script for installing pynauty."""
import platform
import pip
import os
import sys
import tarfile
import shutil
import traceback

from subprocess import check_call

def main(venv=False, cmd_env=False, use_mingw=False):
    """A general main for the installation of pynauty."""

    # Detect OS
    if platform.system() in ['Linux', 'Darwin']:
        check_call(['echo', '"Installing for Unix.."'])
        is_linux = True
        if cmd_env:
            check_call(['echo', "Parameter --cmd_env is ignored.. Installation is happening inside a linus environment"])
        elif use_mingw:
            check_call(['echo', "Parameter --use_mingw is ignored.. Installation is happening inside a linus environment"])

        main_unix(venv)
    elif platform.system() == 'Windows':
        # Warning: For a windows operating system MinGW must be installed.
        # Right now there is no support of pynauty for visual studio
        os.system('echo Installing for windows..')
        if cmd_env or not(use_mingw):
            check_call(['echo', 'Pynauty library is not supported for visual c++ compiler!'])
            check_call(['echo', 'Switching to mingw32..'])

        # Bypass parameter settings
        use_mingw = True
        cmd_env = False

        # Cygwin does not compile either.. Raises error:
        # C:\cygwin64\bin\gcc.exe -mcygwin -mdll -O -Wall -Inauty -Isrc -IC:\Python34\include -IC:\Python34\include -c src/nautywrap.c -o build\temp.win-amd64-3.4\Release\src\nautywrap.o -O4 -fPIC
        # gcc: error: unrecognized command line option '-mcygwin'
        compiler = ' --compiler=mingw32'

        if use_mingw:
            check_call(['echo', 'Importing a mingw32 build library on a python34 system'])
            check_call(['echo', 'the only valid system for intallation of grakel on Windows (because of cvxopt dependencies..)'])
            check_call(['echo', 'allocates a 4 GB memory space on import..'])
            check_call(['echo', 'Even if that is acceptable system-wise, the execution speed will propably decrease..'])
            check_call(['echo', 'This is not a good build..'])
            check_call(['echo', 'For this issue see here:'])
            check_call(['echo', 'https://github.com/ContinuumIO/anaconda-issues/issues/271#issue-58658137'])

        main_windows(venv, cmd_env, use_mingw)
    else:
        sys.stderr.write('Unsupported os for this library')
        sys.exit(1)


def main_windows(venv=False, cmd_env=False, use_mingw=False):
    """Pynauty install on a Windows System"""
    # Add the compiler
    compiler = ' --compiler=mingw32' if use_mingw else ''

    # Running python version
    check_call(['echo', 'python ' + str(platform.python_version())])
    check_call(['echo', 'pip ' + str(pip.__version__)])
    python_executable_address = str(sys.executable)

    # Download library
    check_call(['curl', 'https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz', '-o', 'pynauty-0.6.0.tar.gz'])

    # Decompress
    check_call(['tar', '-xf', 'pynauty-0.6.0.tar.gz'])
    os.remove('pynauty-0.6.0.tar.gz')
    
    # Move inside pynauty
    os.chdir("pynauty-0.6.0")

    # Download C library
    check_call(['curl', 'http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz', '-o', 'nauty26r10.tar.gz'])
    
    # Decompress library
    check_call(['echo', "decompressing: nauty26r10..."])
    check_call(['tar', '-xf', 'nauty26r10.tar.gz'])

    # rename folder
    check_call(['echo', "rename nauty"])
    check_call(['rename', 'nauty26r10', 'nauty'])
    regular_exp = '\"\\\\?\\'+os.path.abspath('nauty')+'\\This_is_nauty26r10.\"'

    # Delete problematic file ending with '.'
    check_call(['echo', "delete config"])
    check_call(['del', regular_exp])

    # build pynauty
    check_call(['echo', 'Making nauty object files ..'])
    check_call(['make', 'nauty-objects'])
    check_call(['echo', 'Object Files made ..'])

    if cmd_env:
        # Set environment variables
        pv = platform.python_version()
        old_pv = os.environ.get('PYTHON_VERSION', '')
        os.environ['PYTHON_VERSION'] = pv

        check_call(['echo', 'Building inside environment'])

        # Execute with the correct environment
        # problem!! nauty.h is designed for gcc compilation
        check_call(['cmd', '/E:ON', '/V:ON', '/C', '..\\ci_scripts\\appveyor\\run_with_env.cmd', python_executable_address, 'setup.py', 'build'])

        # Restore
        if old_pv != '':
            os.environ['PYTHON_VERSION'] = old_pv
    else:
        check_call([python_executable_address, 'setup.py', 'build' + compiler])

    # Install pynauty
    if venv:
        pip.main(['install', '--upgrade', '.'])
    else:
        pip.main(['install', '--user', '--upgrade', '.'])

    # exit directory and delete
    os.chdir("..")
    shutil.rmtree('pynauty-0.6.0')


def main_unix(venv=False):
    """Install on a unix platform"""
    # Running python version
    check_call(['echo', 'python ' + str(platform.python_version())])
    check_call(['echo', 'pip ' + str(pip.__version__)])
    python_executable_address = str(sys.executable)

    # Download library
    check_call(['wget', 'https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz'])

    # Decompress
    check_call(['tar', '-xf', 'pynauty-0.6.0.tar.gz'])
    os.remove('pynauty-0.6.0.tar.gz')

    # Move inside pynauty
    os.chdir("pynauty-0.6.0")

    # Download C library
    check_call(['wget', 'http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz'])

    # Decompress library
    check_call(['echo', "decompressing: nauty26r10..."])
    check_call(['tar', '-xf', 'nauty26r10.tar.gz'])    

    # rename folder
    check_call(['mv', 'nauty26r10', 'nauty'])

    # build pynauty
    check_call(['make', 'nauty-objects'])
    
    # build the executable
    check_call([python_executable_address, 'setup.py', 'build'])

    # Install pynauty
    if venv:
        pip.main(['install', '--upgrade', '.'])
    else:
        pip.main(['install', '--user', '--upgrade', '.'])

    # exit directory and delete
    os.chdir("..")
    shutil.rmtree('pynauty-0.6.0')


if __name__ == "__main__":
    import argparse

    # Create an argument parser for the installer of pynauty
    parser = argparse.ArgumentParser(description='A program to install pynauty')
    parser.add_argument('--venv', help='define if inside a virtual environment', action="store_true")
    if platform.system() == 'Windows':
        meg = parser.add_mutually_exclusive_group()
        meg.add_argument('--cmd_env', help='define if python needs to find the correct environment for building in visual studio [windows]', action="store_true")
        meg.add_argument('--use_mingw', help='define if python needs to use the mingw compiler [windows]', action="store_true")
    args = parser.parse_args()

    # Define parameters useful for installation of pynauty as a consequence of argparse
    argv = [bool(args.venv)]
    if platform.system() == 'Windows':
        argv += [bool(args.cmd_env), bool(args.use_mingw)]

    main(*argv)
