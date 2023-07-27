# ttdilep-delphes

## env setup
```bash
scram p -n delphes CMSSW CMSSW_10_5_0
cd delphes/src
cmsenv
git clone https://github.com/JeewonHeo/ttdilep-delphes.git
mv ttdilep-delphes delphes
scram b -j 20
```

## Delphes 3.4.2
```bash
cmsenv #in the delphes directory
cd ${CMSSW_BASE}/src/delphes/external/
wget http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.4.2.tar.gz
tar -xzvf Delphes-3.4.2.tar.gz
cd Delphes-3.4.2
sed -i s:c++0x:c++17: Makefile
make -j 20
export LD_LIBRARY_PATH=$CMSSW_BASE/src/delphes/external/Delphes-3.4.2:$LD_LIBRARY_PATH
```

## make BuildFile.xml
```bash
cd ${CMSSW_BASE}/src/delphes/analysis/test
./makeBuildFile -d ${CMSSW_BASE}/src/delphes/external/Delphes-3.4.2
```

