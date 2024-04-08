#Uninstall PAPRECA library
cd ../../source/libraries/PAPRECA/
make -f MakeFile clean
rm libKMC.a

cd ../../../Installation/Traditional\ Make/
rm -R ../../TraditionalMakeBuild
