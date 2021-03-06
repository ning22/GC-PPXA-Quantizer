# Copyright 2012 Anna Jezierska
CXX = g++
CC  = g++
# compilator choices

CLINK=-lm 
CXXLINK=-I.
#Makro CLINK: linker options.
GSLLIB = /usr/lib/x86_64-linux-gnu/libgsl.a     # TO BE UPDATED BY THE USER
GSLCBLASLIB = /usr/lib/x86_64-linux-gnu/libgslcblas.a   # TO BE UPDATED BY THE USER

CPPFLAGS=-m64
CFLAGS= -m64 
#Makro CFLAGS: compilator options

CXXSRC =  main.cpp pca.cpp vquan.cpp data_fidelity_term.cpp regularization_term.cpp median_cut.cpp median_cut_vquan.cpp energy.cpp graph_cut_opt.cpp centroid_func.cpp ppxa.cpp graph.cpp maxflow.cpp
COBJ = $(CSRC:.c=.o)
CXXOBJ = $(CXXSRC:.cpp=.o)
OBJ = $(COBJ) $(CXXOBJ)

quantizer: ${CXXSRC} ${CSRC} ${COBJ} ${CXXOBJ}
	$(CXX) ${CPPFLAGS} -o quantizer main.cpp maxflow.o graph.o pca.o vquan.o regularization_term.o data_fidelity_term.o median_cut.o median_cut_vquan.o energy.o graph_cut_opt.o centroid_func.o ppxa.o $(CLINK) $(GSLLIB) $(GSLCBLASLIB)

clean: 
	rm -f $(OBJ) core *.stackdump *.bak

# DO NOT CLEAN
main.o: main.cpp pca.h graph_cut_opt.h centroid_func.h
	$(CXX) -c  $(CPPFLAGS) main.cpp $(GSLLIB) $(GSLCBLASLIB)
pca.o: pca.cpp pca.h
	$(CXX) -c  $(CPPFLAGS) pca.cpp $(GSLLIB) $(GSLCBLASLIB)
vquan.o: vquan.cpp vquan.h median_cut_vquan.h energy.h ppxa.h
	$(CXX) -c  $(CPPFLAGS) vquan.cpp
data_fidelity_term.o: data_fidelity_term.cpp data_fidelity_term.h
	$(CXX) -c  $(CPPFLAGS) data_fidelity_term.cpp
regularization_term.o: regularization_term.cpp regularization_term.h
	$(CXX) -c  $(CPPFLAGS) regularization_term.cpp
median_cut.o: median_cut.cpp median_cut.h
	$(CXX) -c  $(CPPFLAGS) median_cut.cpp	
median_cut_vquan.o: median_cut_vquan.cpp median_cut_vquan.h median_cut.h
	$(CXX) -c  $(CPPFLAGS) median_cut_vquan.cpp	
energy.o: energy.cpp energy.h regularization_term.h data_fidelity_term.h
	$(CXX) -c  $(CPPFLAGS) energy.cpp
maxflow.o: instances.inc graph.h block.h graph.cpp maxflow.cpp 
	$(CXX) -c  $(CPPFLAGS) maxflow.cpp
graph.o: instances.inc graph.h block.h graph.cpp maxflow.cpp 
	$(CXX) -c  $(CPPFLAGS) graph.cpp	
graph_cut_opt.o: graph_cut_opt.cpp graph_cut_opt.h instances.inc energy.h graph.h block.h graph.cpp maxflow.cpp 
	$(CXX) -c  $(CPPFLAGS) graph_cut_opt.cpp
centroid_func.o: centroid_func.cpp centroid_func.h
	$(CXX) -c  $(CPPFLAGS) centroid_func.cpp
ppxa.o: ppxa.cpp ppxa.h
	$(CXX) -c  $(CPPFLAGS) ppxa.cpp $(GSLLIB) $(GSLCBLASLIB)
	


	


