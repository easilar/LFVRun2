rootlibs:=$(shell root-config --libs)
rootflags:=$(shell root-config --cflags)
rootflags2 := -I${ROOTSYS}/include
corlibincl := $(shell correction config --incdir)
#corliblib := /cvmfs/sft.cern.ch/lcg/views/LCG_103/x86_64-centos7-gcc12-opt/lib
corliblib = $(shell a=`correction config --libdir`; echo $${a%%python*lib})
corlibflag := $(shell correction config --cflags --ldflags --rpath)
SOFLAGS = -shared

LD = g++ -m64 -g -Wall

CXXFLAGS = -O0 -g -Wall -fmessage-length=0 $(rootflags) -fpermissive -fPIC -pthread -DSTANDALONE -I. -I$(corlibincl)

OBJDIR=src
SRCDIR=src
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(patsubst %.cpp,%.o,$(SRCS)) $(SRCDIR)/JetMETObjects_dict.o $(SRCDIR)/rootdict.o 

LIBS_EXE = $(rootlibs) -lMathMore -lGenVector -L$(corliblib) -lcorrectionlib
LIBS = $(rootlibs)

TARGET =	nanoaodrdataframe

all:	$(TARGET) libnanoadrdframe.so 

clean:
	rm -f $(OBJS) $(TARGET) libnanoaodrdframe.so $(SRCDIR)/JetMETObjects_dict.C $(SRCDIR)/rootdict.C JetMETObjects_dict_rdict.pcm rootdict_rdict.pcm

$(SRCDIR)/rootdict.C: $(SRCDIR)/NanoAODAnalyzerrdframe.h $(SRCDIR)/TopLFVAnalyzer.h $(SRCDIR)/SkimEvents.h $(SRCDIR)/Linkdef.h
	rm -f $@
	rootcint $@ -I$(corlibincl) -I$(SRCDIR) $^
	rm -f rootdict_rdict.pcm
	ln -s $(SRCDIR)/rootdict_rdict.pcm .

	
libnanoadrdframe.so: $(OBJS)
	$(LD) $(SOFLAGS) $(LIBS_EXE) -D_GLIBCXX_USE_CXX11_ABI=0 -o $@ $^


$(SRCDIR)/JetMETObjects_dict.C: $(SRCDIR)/JetCorrectorParameters.h $(SRCDIR)/SimpleJetCorrector.h $(SRCDIR)/FactorizedJetCorrector.h $(SRCDIR)/JetResolutionObject.h $(SRCDIR)/LinkdefJetmet.h
	rm -f $@
	$(ROOTSYS)/bin/rootcint -f $@ -c $(rootflags2) $^
	ln -s $(SRCDIR)/JetMETObjects_dict_rdict.pcm .

$(SRCDIR)/rootdict.o: $(SRCDIR)/rootdict.C
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(SRCDIR)/JetMETObjects_dict.o: $(SRCDIR)/JetMETObjects_dict.C
	$(CXX) -c -o $@ $(CXXFLAGS) $<

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) -c -o $@ $(CXXFLAGS) $<
	
$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS_EXE)
