# Installation and setup
## Windows only
  Install Windows Subsystem for Linux (WSL) and an Ubuntu version (fx. Ubuntu 20.04 LTS)
  
  https://ubuntu.com/tutorials/install-ubuntu-on-wsl2-on-windows-10#1-overview
  
  (Follow the first 4 steps)
## In WSL/Linux
### Download and configure our git repository
```
git clone https://github.com/JakobScholer/ucloud_jv.git 
```
```
mkdir ucloud_jv/blackbox/zstruct_clones/original/scratch
```
```
mkdir ucloud_jv/blackbox/gsm_clones/original/scratch
```
```
sudo chmod -R 777 ucloud_jv/blackbox/
```
### Install the Intel compiler used by ZStruct and xTB
```
curl -sL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | sudo apt-key add -
```
```
sudo add-apt-repository "deb https://apt.repos.intel.com/oneapi all main"
```
```
sudo apt-get update
```
```
sudo apt-get install intel-hpckit
```
The intel compiler comes with its own python version using conda. We do not want to use conda and disable it by renaming the folder /opt/intel/oneapi/intelpython to disable_intelpython

### Install xTB
Enable the Intel compiler before installing xTB by running
```
source /opt/intel/oneapi/setvars.sh
```
Installing xTB
```
git clone https://github.com/grimme-lab/xtb.git
```
```
cd xtb
```
```
mkdir build
```
```
cd build
```
```
cmake ../ -DCMAKE_BUILD_TYPE=Release
```
```
make
```
```
sudo make install
```
### Install ZStruct
Enable the Intel compiler before installing ZStruct by running
```
source /opt/intel/oneapi/setvars.sh
```
Installing ZStruct
```
git clone https://github.com/JakobScholer/ZStruct.git
```
```
cd ZStruct
```
```
make
```
Copy the created zstruct.exe file into the folder ucloud_jv/blackbox/zstruct_clones/original

### Install Python packages
```
sudo apt install python3-openbabel
```
```
sudo apt install python3-pip
```
```
pip3 install -r requirements.txt
```



# Running program
Individual functions can be called through main.py with arguments, fx:
```
python3 main.py ec
```
```
python3 main.py img_all_stringfiles
```
```
python3 main.py img_stringfile
```
```
python3 main.py make_cut_dag
```
```
python3 main.py smiles_to_reactions_bb
```
```
python3 main.py smiles_to_reactions_nb
```
```
python3 main.py gml
```

# Running tests:
All tests can be run with the command
```
python3 -m unittest discover test
```
Or a single test with the command
```
python3 -m unittest test/filename.py
```
