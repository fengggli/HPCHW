## Compiler, tools and options
CC      = cc
CCFLAGS = -O3 
#CCFLAGS =  -Wall -Wextra  

INCPATH    = -I. -I$(PAPI_ROOT)/include
LIBS       = -lm $(PAPI_LIB)


OBJECTS = summa.o matrix_op.o
TARGET = run_summa



## Implicit rules
.SUFFIXES: .c

.c.o:
	$(CC) -c $(CCFLAGS) $(INCPATH) $<

## Build rules
all:  $(TARGET) 


$(TARGET): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm -f $(OBJECTS)  $(OBJECTS_unroll) $(TARGET_unroll)
	rm -f *~ core

