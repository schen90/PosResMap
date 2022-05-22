#!/bin/bash
#./parallel_AnaData.sh 0 15

Nproc=16  # <----------- max parallel procs

function PushQue {  # push PID into Que
    Que="$Que $1"
    Nrun=$(($Nrun+1))
}

function GenQue {  # update Que
    OldQue=$Que
    Que=""; Nrun=0
    for PID in $OldQue
    do
	if [[ -d /proc/$PID ]]; then
	    PushQue $PID
	fi
    done
}

function ChkQue {  # check Que
    OldQue=$Que
    for PID in $OldQue
    do
	if [[ ! -d /proc/$PID ]]; then
	    GenQue; break
	fi
    done
}

# loop all jobs
for i in `seq $1 $2`
do
    inputfile=rootfiles/G4Sim1/G4SimData$(printf "%04d" $i).root
    outputfile=rootfiles/G4Ana1nc/G4AnaData$(printf "%04d" $i).root

    dbfile=pulsedb/pulseA.root
    #dbfile=pulsedb/SelfCalib/Det0000_PSC.root
    #dbfile=pulsedb/SelfCalib/Det0000_2.0.root
    #dbfile=pulsedb/SelfCalib/Det0000_PSC2.root
    #dbfile=pulsedb/SelfCalib/Det0000_2.0_2.root

    #./macros/AnaData $inputfile $outputfile $dbfile&  # <---------- CMD
    ./macros/GridSearch/build/AnaDataGridSearch $inputfile $outputfile $dbfile&  # <---------- CMD

    PID=$!
    PushQue $PID
    while [[ $Nrun -ge $Nproc ]]
    do
	ChkQue
	sleep 1
    done
done
wait
