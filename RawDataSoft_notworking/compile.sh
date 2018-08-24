#! /bin/bash

g++ -I$ROOTSYS/include -I ../include -fPic -Wall -std=c++11 -m64 -c src/Strip.cc -o src/Strip.o
g++ -Wall -L$ROOTSYS/lib ReadRawData.cc -lCore -lHist src/Strip.o -o main.exe
#g++ -I$ROOTSYS/include -fPic -Wall -c ./src/ReadRawData.cc ./src/Strip.cc -o ReadRawData.exe -std=c++11 `root-config --cflags --glibs --ldflags`
#./src/Event.cc 
