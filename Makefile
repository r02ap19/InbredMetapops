# The directories containing the source files, separated by ':'
CPP_DIRS := src src/* src/*/* 

# A common link flag for all configurations
CPP_LIBPATH := /usr/lib64

#Debug defines
CPP_DEFINES := _DEBUG 

# To make "debug" the default configuration if invoked with just "make":
#
# ifeq ($(CFG),)
# CFG=Debug
# endif

# The source files: regardless of where they reside in the source tree,
CPP_SRC:= $(foreach dir, $(CPP_DIRS), $(wildcard $(dir)/*.cpp))

# Build a Dependency list and an Object list, by replacing the .cpp
# extension to .d for dependency files, and .o for object files.
CPP_DEP = $(patsubst %.cpp, deps.$(CFG)/CPP_%.d, ${CPP_SRC})
CPP_OBJ = $(patsubst %.cpp, objs.$(CFG)/CPP_%.o, ${CPP_SRC})

# Your final binary, executable name
TARGET= your_executable

# What compiler to use for generating dependencies: 
# it will be invoked with -MM -MP
CXXDEP = g++

# What include flags to pass to the compiler
CPP_LIBPATH := $(addprefix -L, $(CPP_LIBPATH))
CPP_INCLUDE := $(addprefix -I, $(CPP_INCLUDE))
CPP_DEFINES := $(addprefix -D, $(CPP_DEFINES))


# Separate compile options per configuration
ifeq ($(CFG),Debug)
CXXFLAGS += -std=c++0x -O0 -g3 -Wall -fpeel-loops ${CPP_DEFINES} ${CPP_INCLUDE} -c
else
CXXFLAGS += -std=c++0x -O3 -Wall -fpeel-loops ${CPP_INCLUDE} -c
endif


all:	inform bin.$(CFG)/${TARGET}

inform:
ifneq ($(CFG),Release)
ifneq ($(CFG),Debug)
	@echo "Invalid configuration "$(CFG)" specified."
	@echo "You must specify a configuration when running make, e.g."
	@echo  "make CFG=Debug"
	@echo  
	@echo  "Possible choices for configuration are 'Release' and 'Debug'"
	@exit 1
endif
endif
	@echo "Configuration "$(CFG)
	@echo "------------------------"

#executable
bin.$(CFG)/${TARGET}: ${CPP_OBJ} | inform
	@mkdir -p $(dir $@)
	$(CXX) -g -o $@ $^ ${CPP_LIBPATH}

#build	
objs.$(CFG)/CPP_%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -MF'$(@:%.o=%.d)' -MT'$(@:%.o=%.d)' -o'$@' '$<'

clean:
	@rm -rf \
	deps.Debug objs.Debug bin.Debug \
	deps.Release objs.Release bin.Release
