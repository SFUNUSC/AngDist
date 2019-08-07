FOR = gfortran -std=legacy -ffixed-line-length-none
LIB = -lmathlib

all:ang_dist ang_dist_cmd

ang_dist: ang_dist.for w.o
	$(FOR) -o ang_dist w.o ang_dist.for $(LIB)
	
ang_dist_cmd: ang_dist_cmd.for w.o
	$(FOR) -o ang_dist_cmd w.o ang_dist_cmd.for $(LIB)

ang_dist_cmd_manual: ang_dist_cmd_manual.for w.o
	$(FOR) -o ang_dist_cmd_manual w.o ang_dist_cmd_manual.for $(LIB)

w.o: w.for
	$(FOR) -c w.for
	
clean:
	rm ang_dist ang_dist_cmd *.o
