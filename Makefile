ROOTCFLAGS   = $(shell root-config --cflags)
ROOTLIBS     = $(shell root-config --libs)
ROOTGLIBS    = $(shell root-config --glibs)
ROOTLDLIB    = $(shell root-config --ldflags)
CXX	     = g++
CXXFLAGS     = -g -Wall $(ROOTCFLAGS)
LDFLAGS      = -g -O2
CINT	     = rootcint
LIBS	     = -lm $(ROOTLIBS) -lRooFit -lRooStats -lTreePlayer -lPhysics $(SYSLIBS) 
BUILD	     = ./build
SRC	     = ./src
INCLUDE      = ./include
LIB	     = ./lib
BIN          = ./bin
RUN	     = ./run
RUN1         = ./run1
RUN2	     = ./run2
RUN3	     = ./run3
PROD0	     = ./prod0
PROD1	     = ./prod1
DICT         = dict
CORE         = core
EXE0         = exe0
EXE1         = exe1

INCFLAGS     = -I$(INCLUDE) 
INCFLAGS    += -I$(TestArea)/DataQuality/GoodRunsLists 
INCFLAGS    += -I$(TestArea)/PhysicsAnalysis/AnalysisCommon/PileupReweighting 
INCFLAGS    += -I$(TestArea)/PhysicsAnalysis/AnalysisCommon/PATCore
INCFLAGS    += -I$(TestArea)/Reconstruction/egamma/egammaAnalysis/egammaAnalysisUtils 
INCFLAGS    += -I$(TestArea)/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection

CXXFLAGS    += $(INCFLAGS)

LOCALLIBS    = -L$(PROD0)/libs -lAnaly -lElectronPhotonFourMomentumCorrection -legammaAnalysisUtils -lGoodRunsLists -lPileupReweighting
EXTLIBS      = -legammaAnalysisUtils -lGoodRunsLists -lPileupReweighting 

CORESRC = $(wildcard $(SRC)/*.cpp)
COREOBJ = $(subst $(SRC)/,$(BUILD)/,$(patsubst %.cpp,%.o,$(CORESRC)) )

EXESRC0     = $(wildcard $(SRC)/$(EXE0)/*.cxx)
EXEOBJ0     = $(subst $(SRC)/$(EXE0)/,$(PROD0)/,$(patsubst %.cxx,%_x,$(EXESRC0)))

EXESRC1     = $(wildcard $(SRC)/$(EXE1)/*.cxx)
EXEOBJ1     = $(subst $(SRC)/$(EXE1)/,$(PROD1)/,$(patsubst %.cxx,%_x,$(EXESRC1)))

default: all

all: obj library cpprod0 cpprod1 exe
all1: obj library cprun1
all2: obj library cprun2
all3: obj library cprun3
prod0: obj library cpprod0
prod1: obj library cpprod1

obj: $(BUILD)/$(DICT)/Dictionary.cpp       $(COREOBJ)

library:
	@echo "in ["$(LIB)"] ---------> Linking library : " $(LIB)/libAnaly.so
	@$(CXX) $(ROOTLDLIB) -shared -o $(LIB)/libAnaly.so  $(COREOBJ) $(BUILD)/$(DICT)/*.o 

exe: $(EXEOBJ0) $(EXEOBJ1) 

cprun:
	@cp $(LIB)/*.so $(RUN)/libs/
	@echo "in ["$(RUN)/libs/"] ---> Copying library : " $(LIB)/libAnaly.so

cprun1:
	@cp $(LIB)/*.so $(RUN1)/libs/
	@echo "in ["$(RUN1)/libs/"] --> Copying library : " $(LIB)/libAnaly.so
cprun2:
	@cp $(LIB)/*.so $(RUN2)/libs/
	@echo "in ["$(RUN2)/libs/"] --> Copying library : " $(LIB)/libAnaly.so
cprun3:
	@cp $(LIB)/*.so $(RUN3)/libs/
	@echo "in ["$(RUN3)/libs/"] --> Copying library : " $(LIB)/libAnaly.so
cpprod0:
	@cp $(LIB)/*.so $(PROD0)/libs/
	@echo "in ["$(PROD0)/libs/"] --> Copying library : " $(LIB)/libAnaly.so
cpprod1:
	@cp $(LIB)/*.so $(PROD1)/libs/
	@echo "in ["$(PROD1)/libs/"] --> Copying library : " $(LIB)/libAnaly.so

$(PROD0)/%_x: $(SRC)/$(EXE0)/%.cxx
	@echo "Building executable : " $@
	@$(CXX) $(CXXFLAGS) -o $@ $? $(LIBS) $(LOCALLIBS)
$(PROD1)/%_x: $(SRC)/$(EXE1)/%.cxx
	@echo "Building executable : " $@
	@$(CXX) $(CXXFLAGS) -o $@ $? $(LIBS) $(LOCALLIBS)

$(BUILD)/$(DICT)/Dictionary.cpp:
	@echo ":::::::::::::: BUILDING ::::::::::::::"
	@echo "in ["$(BUILD)"] -------> Generating global include file and linkdef"
	@macro/buildlinkdef $(INCLUDE) $(BUILD)/$(DICT)/Includes.h $(BUILD)/$(DICT)/Linkdef.h
	@echo "in ["$(BUILD)"] -------> Building ROOT Dictionary."
	@$(CINT) -f $@ -c $(INCFLAGS) $(BUILD)/$(DICT)/Includes.h $(BUILD)/$(DICT)/Linkdef.h
	@$(CXX) $(CXXFLAGS) -I. -fPIC -o $(BUILD)/$(DICT)/Dictionary.o -c $(BUILD)/$(DICT)/Dictionary.cpp




$(BUILD)/%.o: $(SRC)/%.cpp $(INCLUDE)/%.h
	@echo "in ["$(SRC)"] ---------> Building object : " $@
	@$(CXX) $(CXXFLAGS) -fPIC -o $@ -c $<

clean: cleandict cleanobj cleanlib cleanexe

cleanlib:
	@echo "in ["$(LIB)"] -----------> Deleting dynamic library"
	@rm -f $(LIB)/*.so

cleanobj:
	@echo "in ["$(BUILD)"] ---------> Deleting output"
	@rm -f $(BUILD)/*.o
cleandict:
	@echo ":::::::::::::: CLEANING ::::::::::::::"
	@echo "in ["$(BUILD)/$(DICT)"] ----> Deleting CINT dictionnaries"
	@rm -f $(BUILD)/$(DICT)/*.*

cleanexe:
	@echo "in ["$(PROD0)"] ---------> Deleting binaries"
	@rm -f $(PROD0)/*_x
	@echo "in ["$(PROD1)"] ---------> Deleting binaries"
	@rm -f $(PROD1)/*_x
