#! /bin/bash

# build
mkdir -p /app/build
cd /app/build
cmake -G "Unix Makefiles" .. 
cmake --build .

# run
cd /app/bin/Release       
./slide