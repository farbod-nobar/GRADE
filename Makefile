CC=g++
GCCVERS := $(shell gcc -dumpversion)
GCCVERSGT6 := $(shell expr `gcc -dumpversion |cut -f1 -d.` \>= 6 )
ifeq "$(GCCVERSGT6)" "0"
$(info g++ version detected: ${GCCVERS})
$(info g++ version needed >= 6.1.0)
ERR := 121
endif

ifeq "$(ERR)" "121"
$(error GNU version mismatch!)
endif

CFLAGS=-I.

_DEPS = MyFunctions.hpp
DEPS = $(patsubst %,%,$(_DEPS))

_OBJ = GRADE.o MyFunctions.o 
OBJ = $(patsubst %,%,$(_OBJ))


.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

GRADE: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) 

.PHONY: clean

clean:
	rm -f *.o *~ core *~ 
