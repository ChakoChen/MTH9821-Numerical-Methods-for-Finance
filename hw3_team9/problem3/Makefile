# set non-optional compiler flags here
CXXFLAGS += -std=c++1y 

# additional linker flags?
LDFLAGS := -lpthread

# set non-optional preprocessor flags here
# eg. project specific include directories
# CPPFLAGS += 

OBJECTS := Black_Scholes.o Monte_Carlo_Path_Dependent.o Monte_Carlo_Vanilla.o Option.o Random_Number_Generator.o StopWatch.o Test_MC.o 

OUTPUT := MC 

$(OUTPUT) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJECTS) $(LDFLAGS)

Black_Scholes.o : Black_Scholes.hpp Option.hpp Option_Pricer.hpp 

Monte_Carlo_Path_Dependent.o : Monte_Carlo_Path_Dependent.hpp Random_Number_Generator.hpp Option_Pricer.hpp Option.hpp 

Monte_Carlo_Vanilla.o : Monte_Carlo_Vanilla.hpp Random_Number_Generator.hpp Option_Pricer.hpp Option.hpp 

Option.o : Option.hpp 

Random_Number_Generator.o : Random_Number_Generator.hpp 

StopWatch.o : StopWatch.hpp 

Test_MC.o : Option.hpp Black_Scholes.hpp Monte_Carlo_Vanilla.hpp Monte_Carlo_Path_Dependent.hpp Option_Pricer.hpp StopWatch.hpp Random_Number_Generator.hpp 

clean:
	rm $(OUTPUT) $(OBJECTS) *.o
