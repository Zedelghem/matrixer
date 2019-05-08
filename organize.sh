#!/bin/bash
mkdir TO_RUN
for X in nexus/*; do
    # Just the filename, no folder name
    Y=$(echo $X | sed "s/.*\///")
    # Shortest match from the end (.*)
    Z=${Y%.*}
    # Longest match from the end (.*)
    #Z=${Y%%.*}
    
    mkdir TO_RUN/$Z
    cp $X TO_RUN/$Z
    # Change filename to runfile.nexus
    mv TO_RUN/$Z/$Y TO_RUN/$Z/runfile.nexus

# Takes folder path to copy all the files from to every folder in the TO_RUN directory
    if ! [ -z $1 ]; then
        for M in $1/*; do
            cp $M TO_RUN/$Z
        done
    fi

done