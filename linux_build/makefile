CC = g++
DEBUG = -g
CFLAGS =  $(DEBUG) -I../include  -O3
LFLAGS =  -lboost_program_options --static
SRCDIR = ../src/
DEFINES = -D LOGFILE


uai_14: Factor.o graphmodel.o  $(SRCDIR)uai14.cpp
	$(CC) $(CFLAGS)  $(SRCDIR)uai14.cpp graphmodel.o Factor.o $(LFLAGS) -o uai_14



bp_interface: BPinterface.o test_interface.o Factor.o graphmodel.o
	$(CC) $(CFLAGS) $(DEFINES) Factor.o graphmodel.o  BPinterface.o test_interface.o $(LFLAGS)  -o bp_interface

no_opt_bp: BPinterface.o test_interface.o Factor.o graphmodel.o
	$(CC) -g -I../include Factor.o graphmodel.o  BPinterface.o test_interface.o $(LFLAGS)  $(DEFINES) -o bp_interface

test_interface.o: $(SRCDIR)test_interface.cpp
	$(CC) -c -g  -I../include  $(SRCDIR)test_interface.cpp 

BPinterface.o: $(SRCDIR)BPinterface.cpp ../include/BPinterface.h
	$(CC) -c -g $(DEFINES) -I../include  $(SRCDIR)BPinterface.cpp 

graphmodel.o: $(SRCDIR)graphmodel.cpp
	$(CC) -c -g $(DEFINES) -I../include  $(SRCDIR)graphmodel.cpp 

Factor.o: $(SRCDIR)Factor.cpp
	$(CC) -c -g $(DEFINES) -I../include $(SRCDIR)Factor.cpp  



clean:
	\rm *.o *~ uai_14 bp_interface -f 

library_clean:
	\rm *.o bp_interface.class bp_logfile.txt MAR_Sltn.txt libBP.a libBP.so ../test_build/src/libBP.so -f 
