###################################################################################################
# settings
###################################################################################################
HPC_USER      = hpcuser004
HPC_HOME      = /home/$(HPC_USER)

# compiler settings
CC            = gcc
CFLAGS        = -Wall -Ofast -march=native -pipe
LIBS		  = -lm -l OpenCL

# objects
OBJS=heat_diffusion


###################################################################################################
# targets
###################################################################################################
all: clean $(OBJS)

heat_diffusion: heat_diffusion.c
	$(CC) -o $@ $(CFLAGS) $(LIBS) $<

submission: clean heat_diffusion
	tar czf submit.tar.gz heat_diffusion.c Makefile kernel.cl

test: submission
	/home/hpc2019/a4_grading/check_submission.py submit.tar.gz


# helpers
clean:
	rm -f $(OBJS) submit.tar.gz
