#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: ./run.sh <filename.cpp>"
    exit 1
fi

FILE=$1
OUTPUT=${FILE%.*}  # Removes the .cpp extension

g++ -Xpreprocessor -fopenmp "$FILE" -L$(brew --prefix libomp)/lib -lomp -I$(brew --prefix libomp)/include -o "$OUTPUT"
