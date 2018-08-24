
CC = g++
CFLAGS = -I$(IDIR) -I$(ROOTSYS)/include -std=c++11 $(shell root-config --cflags)
LIB = $(shell root-config --libs)
GLIBS  = $(shell root-config --glibs)
IDIR = inc
ODIR = src/obj

_DEPS = ReadRawData.h Event.h Strip.h

DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ReadRawData.o Event.o Strip.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: src/%.cc $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ReadRawData: $(OBJ)
	$(CC) -o main.exe $^ $(CFLAGS) $(LIBS) $(GLIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o main.exe
