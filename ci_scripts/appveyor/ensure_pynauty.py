try:
    import pynauty
    print('Succesful import of pynauty!')
except ImportError as ierr:
    print('Import-Error:', ierr, '..')
    print('Traceback\n' + str(traceback.print_exc()))
except Exception as ex:
	print('Other exception:',ex)