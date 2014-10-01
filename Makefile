FC = ifort
FFLAGS = -Wall -Wextra -r8 -O3
LDFLAGS =
LIBS =
#LIBS = -llapack

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

F90_FILES = $(wildcard src/generate*.f90)
OBJ_FILES = $(patsubst rc/%.f90,obj/%.o,$(F90_FILES))

all: generateSL

generateSL: $(OBJ_FILES)
	$(LINK) -o $@ $^ $(LIBS)

obj/%.o: src/%.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) createSL $(OBJ_FILES) *.mod
