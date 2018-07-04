# set non-optional compiler flags here
CXXFLAGS += -std=c++1y -Wall
#-Wextra -pedantic-errors

# set non-optional preprocessor flags here
# eg. project specific include directories
# CPPFLAGS += 

# find cpp files in subdirectories
SOURCES := $(shell find . -name '*.cpp')

# find headers
HEADERS := $(shell find . -name '*.hpp')

OUTPUT := MC

# Everything depends on the output
all: $(OUTPUT)

# The output depends on sources and headers
$(OUTPUT): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(OUTPUT) $(SOURCES)

clean:
	$(RM) $(OUTPUT) *.o
