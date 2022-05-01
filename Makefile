CXX=clang++
CXXFLAGS=-g -std=c++11 -Wall

all: raytracer2
clean:
	rm -f *.o *.h.gch raytracer2 raytracer2.bmp
test: raytracer2
	./raytracer2
.PHONY: all clean test

raytracer2: raytracer2.o
	$(CXX) $(LDFLAGS) -o $(@) $(^)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $(@) $(<)