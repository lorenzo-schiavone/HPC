#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage: ./run.sh <filename.cpp>"
    exit 1
fi

FILE=$1
OUTPUT=${FILE%.*}  # Removes the .cpp extension

g++ -Xpreprocessor -fopenmp "$FILE" -L/usr/local/opt/libomp/lib -lomp -I/usr/local/opt/libomp/include -o "$OUTPUT"
