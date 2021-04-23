# Install dependencies

```bash
pip3 install --user -r requirements.txt
```

```bash
sudo apt-get install libboost-all-dev libeigen3-dev libhdf5-dev
cd /tmp 
https://github.com/BlueBrain/HighFive.git
cd HighFive
mkdir build && cd build
cmake ..
make && sudo make install
```