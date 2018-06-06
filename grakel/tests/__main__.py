"""The main function for the tests sub-module."""
# Author: Ioannis Siglidis <y.siglidis@gmail.com>
# License: BSD 3 clause

if __name__ == '__main__':
    import os
    import sys
    import warnings

    from subprocess import check_call

    warnings.filterwarnings('ignore', category=UserWarning)

    python_executable_address = str(sys.executable)
    test_dir = str(os.path.dirname(os.path.realpath(__file__)))
    project_dir = str(os.path.realpath(os.path.join(__file__, "../../../")))

    print('Installing the latest "GraKeL"..')
    print('--------------------------------')

    cwd = os.getcwd()
    os.chdir(project_dir)
    try:
        check_call([python_executable_address, project_dir + "/setup.py",
                    "install"])
    finally:
        os.chdir(cwd)

    print('................................................................\n')

    print('Testing Graph..')
    print('---------------')
    check_call([python_executable_address, test_dir + "/test_graph.py",
                "--ignore_warnings", "--verbose"])
    print('................................................................\n')

    print('Testing Kernels..')
    print('-----------------')
    check_call([python_executable_address, test_dir + "/test_kernels.py",
                "--verbose", "--time", "--ignore_warnings", "--all"])
    print('................................................................\n')

    print('Testing Graph Kernels..')
    print('-----------------------')
    check_call([python_executable_address,
                test_dir + "/test_graph_kernels.py",
                "--verbose", "--time", "--ignore_warnings", "--all"])
    print('................................................................')
