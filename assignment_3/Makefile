###################################################################################################
# settings
###################################################################################################
HPC_USER      = hpcuser004
HPC_HOME      = /home/$(HPC_USER)

# compiler settings
CC            = gcc
CFLAGS        = -Wall -Ofast -march=native -pipe
LIBS		  = -lm -lgomp


# objects
OBJS=cell_distances cell_distances.o


###################################################################################################
# targets
###################################################################################################
all: clean $(OBJS)

cell_distances: cell_distances.o
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

cell_distances.o: cell_distances.c
	$(CC) $(CFLAGS) -fopenmp -c $<

submission: clean cell_distances
	tar czf submit.tar.gz cell_distances.c Makefile

test: submission
	/home/hpc2019/a3_grading/check_submission.py submit.tar.gz


# helpers
clean:
	rm -f $(OBJS) submit.tar.gz
