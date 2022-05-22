#!/bin/bash
#./parallel_SimData.sh 300000 1 16

Nproc=8  # <----------- max parallel procs

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
    outputfile=rootfiles/SimData/SimData_run$(printf "%04d" $i).root
    dbfile=pulsedb/pulseA.root
    ./macros/MakeData $1 $dbfile $outputfile &  # <---------- CMD
    PID=$!
    PushQue $PID
    while [[ $Nrun -ge $Nproc ]]
    do
	ChkQue
	sleep 1
    done
done
wait
