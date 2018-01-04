try:
    import pynauty
    print('Succesful import of pynauty!')
except ImportError as ierr:
    print('Import-Error:', ierr, '..')
except Exception as ex:
	print('Other exception:',ex)