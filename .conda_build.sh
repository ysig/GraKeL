set -e
set -x
rm -rf conda_build
mkdir -p conda_build
conda build purge
conda-build . --output-folder conda_build/ --python 2.7
conda-build . --output-folder conda_build/ --python 3.5
conda-build . --output-folder conda_build/ --python 3.6
conda-build . --output-folder conda_build/ --python 3.7
conda-build . --output-folder conda_build/ --python 3.8
conda-build . --output-folder conda_build/ --python 3.9
conda-build . --output-folder conda_build/ --python 3.10
conda-build . --output-folder conda_build/ --python 3.11
conda-build . --output-folder conda_build/ --python 3.12
conda convert -f --platform all conda_build/linux-64/*.tar.bz2 -o conda_build
