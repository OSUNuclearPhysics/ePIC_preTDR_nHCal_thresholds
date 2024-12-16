OBJS = readTreeSimMain.C
EXE = readTreeSimMain

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

INCFLAGS = -I${ROOTSYS}/include -I${ROOTSYS}/include/root/
LDFLAGS = -L${ROOTSYS}/lib -L${ROOTSYS}/lib/root -lpodio -ledm4hep -ledm4eic -lpodioRootIO -lpodioDict -lpodioRootIODict -lfmt

#-L${PODIO}/lib

#CXX = g++ -m32
CXX = g++
FLAGS = -Wall -g $(INCFLAGS) $(LDFLAGS)

COMPILE = $(CXX) $(FLAGS) -c 

all: $(EXE)

#$(EXE): $(OBJS)
#	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS)

$(EXE): $(OBJS)
	$(CXX) -o $(EXE) $(OBJS) $(ROOTFLAGS) $(ROOTLIBS) $(FLAGS)

#%.o: %.C
	$(COMPILE) $<	
	
clean:
	rm $(EXE)

	