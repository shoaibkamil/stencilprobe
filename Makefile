# Makefile for StencilProbe
#
#
CC = gcc
COPTFLAGS = $(PAPI) -O3
CLDFLAGS = $(PAPI)

# the line below defines timers.  if not defined, will attempt to automatically
# detect available timers.  See cycle.h.
# should be set to -DHAVE_PAPI or -DHAVEGETTIMEOFDAY or unset.
#TIMER = -DHAVE_PAPI

probe:	main.c util.c run.h probe_heat.c cycle.h
	$(CC) $(COPTFLAGS) main.c util.c probe_heat.c $(CLDFLAGS) -o probe

blocked_probe:	main.c util.c run.h probe_heat_blocked.c cycle.h
	$(CC) $(COPTFLAGS) main.c util.c probe_heat_blocked.c $(CLDFLAGS) -o probe

clean:
	rm -f *.o probe	
