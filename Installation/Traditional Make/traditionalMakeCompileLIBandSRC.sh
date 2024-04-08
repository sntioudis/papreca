#Install PAPRECA library
cd ../../source/libraries/PAPRECA/
make -f MakeFile clean
make -f MakeFile

cd ../../../Installation/Traditional\ Make/
mkdir ../../TraditionalMakeBuild

make -f MakeFile clean
make -f MakeFile
