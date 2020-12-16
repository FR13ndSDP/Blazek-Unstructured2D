.PHONY : all clean rebuild

MKDIR := mkdir
RM := rm -fr
CC := g++
CXXFLAGS := -O3 -Wall
DIR_DEPS := deps
DIR_EXES := exes
DIR_OBJS := objs

DIRS := $(DIR_DEPS) $(DIR_EXES) $(DIR_OBJS)

EXE := u2d.out
EXE := $(addprefix $(DIR_EXES)/, $(EXE))

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
OBJS := $(addprefix $(DIR_OBJS)/, $(OBJS))
DEPS := $(SRCS:.cpp=.dep)
DEPS := $(addprefix $(DIR_DEPS)/, $(DEPS))

all : $(DIR_OBJS) $(DIR_EXES) $(EXE)

ifeq ("$(MAKECMDGOALS)", "all")
include $(DEPS)
endif

ifeq ("$(MAKECMDGOALS)", "")
include $(DEPS)
endif

$(EXE) : $(OBJS)
	$(CC) $^ $(CXXFLAGS) -o $@ 
	@echo "Success! Target => $@"

$(DIR_OBJS)/%.o : %.cpp
	$(CC) $(CXXFLAGS) -o $@ -c $(filter %.cpp, $^)
#   $(CC) -o $@ -c $(filter %.c, $^)

$(DIRS) :
	$(MKDIR) $@

ifeq ("$(wildcard $(DIR_DEPS))", "")
$(DIR_DEPS)/%.dep : $(DIR_DEPS) %.cpp
else
$(DIR_DEPS)/%.dep : %.cpp
endif
	@echo "Creating $@ ..."
	@set -e; \
	$(CC) -MM -E $(filter %.cpp, $^) | sed 's,\(.*\)\.o[ :]*,objs/\1.o $@ : ,g' > $@

clean :
	$(RM) $(DIRS)

rebuild :
	@$(MAKE) clean
	@$(MAKE) all
