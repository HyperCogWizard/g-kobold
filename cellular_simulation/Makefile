# Cellular Simulation Module Makefile

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra -I.
LDFLAGS = 

# Source files
CORE_SOURCES = core/atomspace.cpp core/psystem.cpp
MEMBRANE_SOURCES = 
ALL_SOURCES = $(CORE_SOURCES) $(MEMBRANE_SOURCES)

# Object files
CORE_OBJECTS = $(CORE_SOURCES:.cpp=.o)
MEMBRANE_OBJECTS = $(MEMBRANE_SOURCES:.cpp=.o)
ALL_OBJECTS = $(ALL_SOURCES:.cpp=.o)

# Test executable
TEST_TARGET = test_cellular
TEST_SOURCES = test_cellular.cpp

.PHONY: all clean test

all: $(TEST_TARGET)

$(TEST_TARGET): $(ALL_OBJECTS) $(TEST_SOURCES)
	$(CXX) $(CXXFLAGS) -o $@ $(ALL_OBJECTS) $(TEST_SOURCES) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test: $(TEST_TARGET)
	./$(TEST_TARGET)

clean:
	rm -f $(ALL_OBJECTS) $(TEST_TARGET)

# Dependencies
core/atomspace.o: core/atomspace.cpp core/atomspace.h
core/psystem.o: core/psystem.cpp core/psystem.h core/atomspace.h

# Show what we're building
info:
	@echo "Core sources: $(CORE_SOURCES)"
	@echo "All objects: $(ALL_OBJECTS)"
	@echo "Test target: $(TEST_TARGET)"