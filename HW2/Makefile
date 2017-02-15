## Compiler, tools and options
CC      = cc
CCFLAGS = -O3 
#CCFLAGS =  -Wall -Wextra  

## PAPI Setup
PAPI_ROOT  = /opt/cray/papi/5.3.0
PAPI_LIB   = $(PAPI_ROOT)/lib/libpapi.a

INCPATH    = -I. -I$(PAPI_ROOT)/include
LIBS       = -lm $(PAPI_LIB)


OBJECTS_unroll = run_unroll.o  papi_timer.o
TARGET_unroll  = run_unroll



## Implicit rules
.SUFFIXES: .c

.c.o:
	$(CC) -c $(CCFLAGS) $(INCPATH) $<

## Build rules
all:  $(TARGET_unroll) 


$(TARGET_unroll): $(OBJECTS_unroll)
	$(CC) -o $@ $(OBJECTS_unroll) $(LIBS)

clean:
	rm -f $(OBJECTS)  $(OBJECTS_unroll) $(TARGET_unroll)
	rm -f *~ core

