#**********************************************************************
# Author:		Matej Pavin
# Email			mpavin@triumf.ca
# Date:			02/16/2018
# Experiment:	Emphatic (T1396)
#
# Compile. Or do not compile. There is no try. - Yoda
#********************************************************************** 
INCDIR= $(shell pwd)/include
SRCDIR= $(shell pwd)/src
OBJDIR= $(shell pwd)/obj
BINDIR= $(shell pwd)/bin

VPATH = $(SRCDIR)

CXX=g++
CFLAGS=-c -O3 -fopenmp -g -Wall `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter -lMinuit -ltbb
#LDFLAGS=`root-config --glibs` -lHistPainter -lMinuit -lgomp


TARGET=waveform_analysis.cc

EXECUTABLE=$(TARGET:%.cc=$(BINDIR)/%.app)

FILES= $(wildcard $(SRCDIR)/*.cc)
SOURCES=$(FILES)

OBJECTS = $(FILES:$(SRCDIR)/%.cc=${OBJDIR}/%.o)

OBJ=$(TARGET:%.cc=${OBJDIR}/%.o) $(OBJECTS)

all: MESSAGE $(TARGET) $(SOURCES) $(EXECUTABLE)

MESSAGE:
	@echo '**********************************************************************'
	@echo '* Author:     Matej Pavin                                            *'
	@echo '* Email       mpavin@triumf.ca                                       *'
	@echo '* Date:       2022/07/14			                            *'
	@echo '* Water Cherenkov Test Experiment (WCTE)                             *'
	@echo '**********************************************************************'	

$(EXECUTABLE): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CFLAGS) $< -o $@

print-%  : ; @echo $* = $($*)

clean:
	- $(RM) $(BINDIR)/* $(OBJDIR)/*
