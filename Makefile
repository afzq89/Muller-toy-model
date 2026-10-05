OBJECTS = toy.o
CCOMPLR = g++ -O2
CCFLAGS = -g -std=gnu++0x

toyUS: toyUS.cc
	$(CCOMPLR) $(CCFLAGS) toyUS.cc -o $@
