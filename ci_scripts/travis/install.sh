# Install file used by cibuildwheel
pip install --upgrade pip
pip install --upgrade setuptools
pip install -r requirements.txt
pip install "cvxopt>=1.2.0"

# Install coverage if needed
if [[ "$COVERAGE" == "true" ]]; then
    pip install networkx
    pip install pandas
    pip install coverage
fi

# Print versions
python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import cython; print('cython %s' % cython.__version__)"
python -c "import cvxopt; print('cvxopt %s' % cvxopt.__version__)"

