#!/usr/bin/env sh

# Requirements
# sudo apt-get install python-dev
# sudo apt-get install python-dev
# sudo apt-get install pip3
# pip3 install setuptools
wget https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz
tar -xf pynauty-0.6.0.tar.gz
rm -rf pynauty-0.6.0.tar.gz
cd pynauty-0.6.0
wget http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz
tar -xf nauty26r10.tar.gz
ln -s nauty26r10 nauty
make pynauty
make tests
sudo make user-ins
cd ..
rm -rf pynauty-0.6.0
