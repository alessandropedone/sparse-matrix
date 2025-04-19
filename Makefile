# Compiler flags
CXX      ?= g++
CXXFLAGS ?= -std=c++20 -Wall -O3
CPPFLAGS ?= -I include  # Include flags

# Linker flags
LDFLAGS ?=
LDLIBS  ?= -ltbb

# Variables
EXEC    = main
SRC_DIR = src
SRCS 	= $(shell find $(SRC_DIR) -name '*.cpp')
OBJS    = $(SRCS:.cpp=.o)
HEADERS = $(shell find include -maxdepth 1 -name '*.hpp')

# Default target
all: $(EXEC)

# Link object files to create executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $@

# Compile source files
%.o: %.cpp $(HEADERS)
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

# Remove all object files
clean:
	$(RM) $(OBJS)

# Remove all generated files
distclean: clean
	$(RM) $(EXEC)
	$(RM) $(SRC_DIR)/*~