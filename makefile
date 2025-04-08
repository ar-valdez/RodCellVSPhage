PROGRAM := Main

COMPILER := g++ -std=c++17

#	Directories of source code and object code
SRCDIR    = .
OBJDIR    = obj
TARGETDIR = run_sim

#	$< stands for input or dependency
#	$@ stands for output or target
#	$^ stands for inputs

SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
LIBS     := -fopenmp -O3

#	Linking
$(TARGETDIR)/$(PROGRAM): $(OBJECTS)
	$(COMPILER) $^ -o $@ $(LIBS)

#	Compiling
$(OBJECTS): $(OBJDIR)/%.o: %.cpp
	$(COMPILER) -c $< -o $@ $(LIBS)

clean:
	rm -f $(OBJDIR)/*.o
