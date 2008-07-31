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
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat.c $(CLDFLAGS) -o probe

circqueue_probe:	main.c util.c run.h probe_heat_circqueue.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES -DCIRCULARQUEUEPROBE main.c util.c probe_heat_circqueue.c $(CLDFLAGS) -o probe

timeskew_probe:	main.c util.c run.h probe_heat_timeskew.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_timeskew.c $(CLDFLAGS) -o probe

oblivious_probe:	main.c util.c run.h probe_heat_oblivious.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_oblivious.c $(CLDFLAGS) -o probe

blocked_probe:	main.c util.c probe_heat_blocked.c cycle.h
	$(CC) $(COPTFLAGS) $(TIMER) -DRANDOMVALUES main.c util.c probe_heat_blocked.c $(CLDFLAGS) -o probe

test:	main.c util.c run.h probe_heat.c cycle.h  probe_heat_blocked.c probe_heat_oblivious.c probe_heat_timeskew.c probe_heat_circqueue.c
	$(CC) $(COPTFLAGS) -DSTENCILTEST main.test.c util.c probe_heat.c probe_heat_blocked.c probe_heat_oblivious.c probe_heat_timeskew.c probe_heat_circqueue.c $(CLDFLAGS) -o probe

clean:
	rm -f *.o probe	
