CC=g++
CFLAGS=-fPIC
INC=-I$(HOME)/src/OPCODE/lib/HP2C_DYCORE/src/stencil_framework -I$(HOME)/src/OPCODE/lib/HP2C_DYCORE/src/shared_infrastructure -I$(HOME)/src/OPCODE/lib/HP2C_DYCORE/src/shared_definitions -I$(HOME)/src/OPCODE/lib/HP2C_DYCORE/src/hp2c_dycore -I./include/python3.4m
LIB=libHP2CDycore.a libSharedInfrastructure.a
OBJ = Coriolis.o stella_wrap.o
BIN = _stella.so

%.o: %.cpp
		$(CC) -c -o $@ $< $(CFLAGS) $(INC)

$(BIN): swig $(OBJ)
		$(CC) -shared -o $@ $(OBJ) $(CFLAGS) $(LIB)

swig: 
		swig -c++ -python -py3 stella.i
		mv stella_wrap.cxx stella_wrap.cpp

clean:
	rm -rf $(DEPS) $(OBJ) $(BIN)
