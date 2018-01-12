if __name__=='__main__':
    import os, sys
    from subprocess import call
    import warnings

    warnings.filterwarnings('ignore', category=UserWarning)

    python_executable_address = str(sys.executable)
    test_dir = str(os.path.dirname(os.path.realpath(__file__)))
    project_dir = str(os.path.abspath(os.path.join(__file__ ,"../../../")))


    print('Installing the latest "GraKeL"..')
    print('------------------------------')
    call([python_executable_address, project_dir + "/setup.py", "install"])
    print('................................................................\n')

    print('Testing Graph..')
    print('---------------')
    call([python_executable_address, test_dir + "/test_graph.py", "--ignore_warnings", "--verbose"])
    print('................................................................\n')

    print('Testing Kernels..')
    print('-----------------')
    call([python_executable_address, test_dir + "/test_kernels.py", "--verbose", "--ignore_warnings", "--all"])
    print('................................................................\n')

    print('Testing Graph Kernels..')
    print('-----------------------')
    call([python_executable_address, test_dir + "/test_graph_kernels.py", "--verbose", "--ignore_warnings", "--all"])
    print('................................................................')