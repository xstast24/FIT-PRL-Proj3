#!/bin/bash
# Merge Splitting Sort (PRL project 2)
# Author: Filip Stastny (xstast24)

if [ $# != 1 ]; then
    echo "Ocekava 1 parametr: retezec uzlu bin. stromu"
    exit 1
else
    retezec=$1;
    delka_retezce=${#retezec};
    let proc_count=2*$delka_retezce-2;
    if [ $proc_count == 0 ]; then # in case of a single item, run only one processor
        proc_count=1;
    fi
fi

mpic++ --prefix /usr/local/share/OpenMPI -o proj3 proj3.cpp
mpirun --prefix /usr/local/share/OpenMPI -np $proc_count proj3 $retezec
# lcoal values
# mpic++ -o proj3 proj3.cpp
# mpirun -np $proc_count proj3 $retezec