## Compiler, tools and options
CC      = cc
CCFLAGS = -O3 -Wall -Wextra  -g

## PAPI Setup
PAPI_ROOT  = /opt/cray/papi/5.3.0
PAPI_LIB   = $(PAPI_ROOT)/lib/libpapi.a

INCPATH    = -I. -I$(PAPI_ROOT)/include
LIBS       = -lm $(PAPI_LIB)

## Files
OBJECTS_block = run_block.o  papi_timer.o
TARGET_block  = run_block

OBJECTS_unroll = run_unroll.o  papi_timer.o
TARGET_unroll  = run_unroll



## Implicit rules
.SUFFIXES: .c

.c.o:
	$(CC) -c $(CCFLAGS) $(INCPATH) $<

## Build rules
all:  $(TARGET_unroll) $(TARGET_block)

$(TARGET_block): $(OBJECTS_block)
	$(CC) -o $@ $(OBJECTS_block) $(LIBS)

$(TARGET_unroll): $(OBJECTS_unroll)
	$(CC) -o $@ $(OBJECTS_unroll) $(LIBS)

clean:
	rm -f $(OBJECTS)  $(OBJECTS_block) $(TARGET_block) $(OBJECTS_unroll) $(TARGET_unroll)
	rm -f *~ core

