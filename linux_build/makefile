CC = g++
DEBUG = -g
CFLAGS =  $(DEBUG) -I../include  -O3
LFLAGS =  -lboost_program_options --static
SRCDIR = ../src/

uai_14: $(SRCDIR)Factor.cpp $(SRCDIR)factorgraph.cpp $(SRCDIR)graphmodel.cpp  $(SRCDIR)uai14.cpp
	$(CC) $(CFLAGS) $(SRCDIR)Factor.cpp $(SRCDIR)factorgraph.cpp $(SRCDIR)graphmodel.cpp  $(SRCDIR)uai14.cpp $(LFLAGS) -o uai_14

bp_interface: $(SRCDIR)BPinterface.cpp $(SRCDIR)test_interface.cpp
	$(CC) $(CFLAGS) $(SRCDIR)Factor.cpp $(SRCDIR)factorgraph.cpp $(SRCDIR)graphmodel.cpp  $(SRCDIR)BPinterface.cpp $(SRCDIR)test_interface.cpp $(LFLAGS)  -o bp_interface

no_opt_bp: $(SRCDIR)BPinterface.cpp $(SRCDIR)test_interface.cpp
	$(CC) -g -I../include  $(SRCDIR)Factor.cpp $(SRCDIR)factorgraph.cpp $(SRCDIR)graphmodel.cpp  $(SRCDIR)BPinterface.cpp $(SRCDIR)test_interface.cpp $(LFLAGS)  -o bp_interface

clean:
	\rm *.o *~ uai_14 bp_interface -f

