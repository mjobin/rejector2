CXXFLAGS=-O2 -fpermissive -pthread


rejector2source := $(wildcard src/rejector2/*.cpp)
rejector2objects := $(patsubst %.cpp,%.o,$(wildcard src/rejector2/*.cpp))


rej2runsource := $(wildcard src/rej2run/*.cpp)
rej2runobjects := $(patsubst %.cpp,%.o,$(wildcard src/rej2run/*.cpp))

all : rejector2 rej2run


	
rejector2 : $(rejector2objects)
	$(CXX) $(CXXFLAGS) -o rejector2 $(rejector2objects)



rej2run : $(rej2runobjects)
	$(CXX) $(CXXFLAGS) -o rej2run $(rej2runobjects)
