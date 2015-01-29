CC = g++
DEBUG = -g
CFLAGS =  $(DEBUG) -I../include	 -O3
SRCDIR = ../src/


bp_interface: BPinterface.o test_interface.o Factor.o graphmodel.o
	$(CC) $(CFLAGS) -D WINDOWS	Factor.o graphmodel.o  BPinterface.o test_interface.o  -o bp_interface

bp_library: BPinterface.o graphmodel.o Factor.o
	g++ -shared -fPIC -D WINDOWS -o LBP.dll BPinterface.o graphmodel.o Factor.o

april_test: APPRIL_lib_test.cpp
	g++ -g APPRIL_lib_test.cpp  -I../include -L./LBP.dll -o april_test
	
	
test_interface.o: $(SRCDIR)test_interface.cpp
	$(CC) -c -g -D WINDOWS	 -I../include  $(SRCDIR)test_interface.cpp 

BPinterface.o: $(SRCDIR)BPinterface.cpp ../include/BPinterface.h
	$(CC) -c -g  -D WINDOWS	 -I../include  $(SRCDIR)BPinterface.cpp 

graphmodel.o: $(SRCDIR)graphmodel.cpp
	$(CC) -c -g -D WINDOWS	-I../include  $(SRCDIR)graphmodel.cpp 

Factor.o: $(SRCDIR)Factor.cpp
	$(CC) -c -g  -D WINDOWS	 -I../include $(SRCDIR)Factor.cpp  



clean:
	\rm *.o *~ uai_14 bp_interface -f LBP.dll MAR_Sltn.txt bp_logfile.txt
	

