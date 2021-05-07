###################################################################################################
# settings
###################################################################################################
HPC_USER      = hpcuser004
HPC_HOME      = /home/$(HPC_USER)

# compiler settings
CC            = mpicc
CFLAGS        = -Wall -Ofast -march=native -pipe
LIBS		  = -lm

# objects
OBJS=heat_diffusion


###################################################################################################
# targets
###################################################################################################
all: clean $(OBJS)

heat_diffusion: heat_diffusion.c
	$(CC) -o $@ $(CFLAGS) $(LIBS) $<

submission: clean heat_diffusion
	tar czf submit.tar.gz heat_diffusion.c Makefile

test: submission
	/home/hpc2019/a5_grading/check_submission.py submit.tar.gz


# helpers
clean:
	@ rm -f $(OBJS) submit.tar.gz
