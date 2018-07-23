PROJ_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ifndef OCCA_DIR
  include $(PROJ_DIR)/../../../scripts/Makefile
else
  include ${OCCA_DIR}/scripts/Makefile
endif


paths += -I./include 
flags += -DOCCA_GL_ENABLED=1 -Dp_N=$(N) -g 

ifeq ($(OS),OSX)
	links += -framework OpenGL -framework GLUT
endif

ifeq ($(OS),LINUX)
	links +=
# -lGLU -lglut
endif


#---[ COMPILATION ]-------------------------------
headers = $(wildcard $(incPath)/*.hpp) $(wildcard $(incPath)/*.tpp)
sources = $(wildcard $(srcPath)/*.cpp)

objects  = $(subst $(srcPath)/,$(objPath)/,$(sources:.cpp=.o))

.PHONY: clean

all: main2d main3d

main2d: $(objects) $(headers) main2d.cpp
	$(compiler) $(compilerFlags) -o main2d $(flags) $(objects) main2d.cpp -L${OCCA_DIR}/lib $(paths) $(links)

main3d: $(objects) $(headers) main3d.cpp
	$(compiler) $(compilerFlags) -o main3d $(flags) $(objects) main3d.cpp -L${OCCA_DIR}/lib $(paths) $(links)

$(objPath)/%.o:$(srcPath)/%.cpp $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.tpp)))
	$(compiler) $(compilerFlags) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f $(objPath)/*.o;
	rm -f main;
	rm -f main2d;
	rm -f main3d;
	rm -rf main2d.dSYM/;
	rm -rf main3d.dSYM/;

#=================================================
