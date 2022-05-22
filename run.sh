#!/bin/bash
#./run.sh 1

#obj=/home/userfs/s/sc2950/server/agata/trunk/build-gcc9.3.0/Agata
obj=/home/sc2950/agata/trunk/build-gcc9.3.0/Agata
mac=G4mac/sim1gamma.mac

#($obj -seed -Path ./trunk/ -b $mac -run $1)
($obj -Path ./trunk/ -b $mac -run $1)

./macros/MakeDataG4 trunk/GammaEvents.$(printf "%04d" $1) rootfiles/G4Sim/G4SimData$(printf "%04d" $1).root
