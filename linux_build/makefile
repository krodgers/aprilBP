CC = g++
DEBUG = -g
CFLAGS =  $(DEBUG) -I../include  -O3
LFLAGS =  -lboost_program_options --static
SRCDIR = ../src/


uai_14: factor.o graphmodel.o  $(SRCDIR)uai14.cpp
	$(CC) $(CFLAGS)  $(SRCDIR)uai14.cpp graphmodel.o factor.o $(LFLAGS) -o uai_14


take_two: BpInterfaceTakeTwo.o test_interface.o factor.o graphmodel.o
	$(CC) $(CFLAGS) factor.o graphmodel.o  BPinterface.o test_interface.o $(LFLAGS)  -o bp_interface

bp_interface: BPinterface.o test_interface.o factor.o graphmodel.o
	$(CC) $(CFLAGS) factor.o graphmodel.o  BPinterface.o test_interface.o $(LFLAGS)  -o bp_interface

no_opt_bp: BPinterface.o test_interface.o factor.o graphmodel.o
	$(CC) -g -I../include factor.o graphmodel.o  BPinterface.o test_interface.o $(LFLAGS)  -o bp_interface

test_interface.o: $(SRCDIR)test_interface.cpp
	$(CC) -c -g  -I../include  $(SRCDIR)test_interface.cpp 

BPinterface.o: $(SRCDIR)BPinterface.cpp ../include/BPinterface.h
	$(CC) -c -g -I../include  $(SRCDIR)BPinterface.cpp 

BpInterfaceTakeTwo.o: $(SRCDIR)BpInterfaceTakeTwo.cpp ../include/BpInterfaceTakeTwo.h
	$(CC) -c -g -I../include  $(SRCDIR)BpInterfaceTakeTwo.cpp 

graphmodel.o: $(SRCDIR)graphmodel.cpp
	$(CC) -c -g -I../include  $(SRCDIR)graphmodel.cpp 

factor.o: $(SRCDIR)Factor.cpp
	$(CC) -c -g -I../include $(SRCDIR)Factor.cpp  



clean:
	\rm *.o *~ uai_14 bp_interface -f

