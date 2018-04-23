# Deactivate the travis-provided virtual environment and setup a
# conda-based environment instead
deactivate

# Use the miniconda installer for faster download / install of conda
# itself
pushd .
cd
mkdir -p download
cd download
echo "Cached in $HOME/download :"
ls -l
echo

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
   wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
else
   wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
fi

chmod +x miniconda.sh && ./miniconda.sh -b
cd ..

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
   export PATH=/home/travis/miniconda2/bin:$PATH;
else
   export PATH=/home/travis/miniconda3/bin:$PATH;
fi

conda update --yes conda
popd

# Configure the conda environment and put it in the path using the
# provided versions
conda create -n testenv --yes python=$TRAVIS_PYTHON_VERSION pip nose \
      numpy=$NUMPY_VERSION scipy=$SCIPY_VERSION cython=$CYTHON_VERSION
      
source activate testenv

if [[ "$COVERAGE" == "true" ]]; then
    pip install coverage coveralls
fi

pip install 'cvxopt>=1.2.0'
pip install flake8
# If setuptools are not installed -- upgraded, travis crashes
pip install --upgrade setuptools

python --version
python -c "import numpy; print('numpy %s' % numpy.__version__)"
python -c "import scipy; print('scipy %s' % scipy.__version__)"
python -c "import cvxopt; print('cvxopt %s' % cvxopt.__version__)"
python setup.py develop
