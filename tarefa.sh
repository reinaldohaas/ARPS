   sudo apt update
   sudo apt upgrade
   sudo apt install gfortran libnetcdff-dev libhdf4-dev
   sudo apt install libhdf4-dev
   sudo apt install libhdf5-dev
   sudo apt install openmpich
   sudo apt install libmpich-dev
   wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB | gpg --dearmor | sudo tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
    echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | sudo tee /etc/apt/sources.list.d/oneAPI.list
   sudo apt update
   sudo apt install intel-oneapi-hpc-toolkit
