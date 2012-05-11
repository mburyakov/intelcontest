# sample Makefile.
# It compiles every .cpp files in the src/ directory to object files in obj/ directory, and build the ./run executable.
# It automatically handles include dependencies.

# You can modify it as you want (add libraries, compiler flags...).
# But if you want it to properly run on our benchmark service, don't rename or delete variables.

# using icc :
#COMPILER ?= $(ICC_PATH)icpc
# using gcc :
COMPILER ?= $(GCC_PATH)g++

# using icc :
#FLAGS ?= -std=c++0x -U__GXX_EXPERIMENTAL_COMPILER0X__ -xHOST -fast -w1 $(ICC_SUPPFLAGS)
# using gcc :
FLAGS ?= -std=c++0x -O3 -Wall $(GCC_SUPPFLAGS)

LDFLAGS ?= -g
LDLIBS =  -lipps_l -lippcore_l
INCLUDEPATHS = -I/home/mburyakov/bin/intel/composer_xe_2011_sp1.10.319/ipp/include/ -L/home/mburyakov/bin/intel/composer_xe_2011_sp1.10.319/ipp/lib/intel64
#example if using Intel� Threading Building Blocks :
#LDLIBS = -ltbb -ltbbmalloc

EXECUTABLE = run
ARGS = 24 16 refseq.txt input.txt

TEAM_ID = 19e8fda99f85c066c1cc5e7dcee5f86f # put your 32chars team id here and you will be able to submit your program from command line using "make submit"

SRCS=$(wildcard src/*.cpp)
OBJS=$(SRCS:src/%.cpp=obj/%.o)

all: release

release: $(OBJS)
	$(COMPILER) $(LDFLAGS)  -o $(EXECUTABLE) $(OBJS) $(INCLUDEPATHS) $(LDLIBS)

obj/%.o: src/%.cpp
	$(COMPILER) $(FLAGS) -o $@ -c $< $(INCLUDEPATHS)

run: all
	./$(EXECUTABLE) $(ARGS)


zip: dist-clean
ifdef TEAM_ID
	zip $(strip $(TEAM_ID)).zip -r ./
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif

submit: zip
ifdef TEAM_ID
	curl -F "file=@$(strip $(TEAM_ID)).zip" -L http://intel-software-academic-program.com/ayc-upload/upload.php
else
	@echo "you need to put your TEAM_ID in the Makefile"
endif

clean:
	rm -f obj/*

dist-clean: clean
	rm -f $(EXECUTABLE) *~ .depend *.zip
	
#automatically handle include dependencies
depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	@$(foreach SRC, $(SRCS), $(COMPILER) $(FLAGS) -MT $(SRC:src/%.cpp=obj/%.o) -MM $(SRC) >> .depend;)

include .depend	
