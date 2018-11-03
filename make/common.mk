# Copyright (c) 2017-2019 Wendy S. Alfieri
#
# WARNING: THIS DOCUMENT CONTAINS TRADE SECRET INFORMATION OF WENDY S. ALFIERI.  
# UNAUTHORIZED DISCLOSURE IS STRICTLY PROHIBITED AND MAY RESULT IN SERIOUS LEGAL CONSEQUENCES.
#
# O/S Dependencies and Common Rules for Makefiles - the disgusting stuff
#
OS = $(shell uname)

############################
# DEFAULTS - gcc-based toolset
############################

CC      = gcc
LD      = gcc
CFLAGS  = -std=c++17 -Wextra -Wstrict-aliasing -pedantic -fmax-errors=10 -Werror -Wunreachable-code -Wcast-align -Wcast-qual
CFLAGS += -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-include-dirs -Wnoexcept 
CFLAGS += -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-promo -Wstrict-null-sentinel 
#CFLAGS += -Wstrict-overflow=5
CFLAGS += -Wswitch-default -Wundef 
CFLAGS += -fno-strict-aliasing
CFLAGS += -fdiagnostics-color=never 
CFLAGS += -mavx2 -mfma -march=haswell -g -O3
CFLAGS += -I. -I../lib -I/home/tools/lapack-3.7.1/LAPACKE/include
LFLAGS  = -g -lm -lz -lstdc++ -lpthread -rdynamic -fdiagnostics-color=always 
DRUN    = gdb -tui


############################
# LINUX OVERRIDES
############################
ifeq ($(OS), Linux)

GLUT_DIR = /home/utils/freeglut-2.8.1
CFLAGS = -Wall -Werror -pedantic -Wno-long-long -Wno-deprecated -O3 -g -DEMULATE_BUFFERS -I../base -I${GLUT_DIR}/include
LFLAGS = -g -lm -lz -lstdc++ -lGL -lglut -lGLU -L${GLUT_DIR}/lib

############################
# MACOS OVERRIDES
############################
else
ifeq ($(OS), Darwin)

CC = clang
LD = clang -w
CFLAGS += -Wno-unknown-warning-option -Wno-unused-command-line-argument -Wno-unused-parameter -Wno-shift-count-negative -DEMULATE_BUFFERS -DGLUT_ONLY  -I../base -I/usr/include/malloc/ -DNO_MALLOC_H -DUSE_POSIX_MEMALIGN -Wno-deprecated-declarations
LFLAGS = -g -lm -lz -lstdc++ -framework OpenGl -framework GLUT -framework CoreFoundation -framework IOKit -framework Carbon -framework CoreGraphics -framework Cocoa
DRUN   = lldb -s .gdbinit

############################
# CYGWIN OVERRIDES
############################
else
ifneq (,$(findstring CYGWIN, $(OS)))

CFLAGS = -Wall -Werror -pedantic -O3 -g -DEMULATE_BUFFERS -I../base
LFLAGS = -g -lm -lz -lstdc++ -lGL -lglu32 -lglut

else
$(error Unknown O/S: $(OS))

endif
endif
endif

############################
# COMMON RULES
############################

CFLAGS += $(EXTRA_CFLAGS)
LFLAGS += $(EXTRA_LFLAGS)

RUNS = $(PROGS:.exe=.run)

all: $(OBJS) $(PROGS)

check: $(RUNS)

.PRECIOUS: %.o

# compile
%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $<

# link
%.exe: %.o $(DEPS) $(OBJS) $(LIBS)
	$(LD) $(OBJS) $*.o $(LIBS) $(LFLAGS) -o $@ 

# normal run
%.run: %.exe
	./$<

# debugger
%.drun: %.exe
	$(DRUN) $*.exe

regress:
	@echo "Building everything..."; \
	make; \
        echo ""; \
	for f in $(REGRESSION); do \
            echo -n "Running $${f}..."; \
	    make $$f.run > /dev/null; \
            rc=$$?; if [ $$rc != 0 ]; then exit $$rc; fi; \
	    echo " PASSED"; \
	done; \
        echo ""; \
	echo "ALL PASSED"

clean:
	rm -fr *.o *.s *.exe .libs *.gz *.out ._* .libs $(PROGS) *.stackdump *.mp4 *.cxx *.vcmc *.vmc

