# Configure the conda environment and put it in the path using the
# provided versions
pip install nose "numpy==$NUMPY_VERSION" "scipy==$SCIPY_VERSION" "cython==$CYTHON_VERSION" future six 'cvxopt>=1.2.0' flake8 --upgrade

if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage coveralls
fi

# If setuptools are not installed -- upgraded, travis crashes
pip install --upgrade setuptools

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import cython; print('cython %s' % cython.__version__)"
python -c "import cvxopt; print('cvxopt %s' % cvxopt.__version__)"
python setup.py develop
