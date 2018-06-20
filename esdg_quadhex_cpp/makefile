PROJ_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ifndef OCCA_DIR
  include $(PROJ_DIR)/../../../scripts/Makefile
else
  include ${OCCA_DIR}/scripts/Makefile
endif


paths += -I./include 
flags += -DOCCA_GL_ENABLED=1 -Dp_N=$(N) -g -L${OCCA_DIR}/lib

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

# curved elements with WADG (can also add heterogeneous)
main: $(objects) $(headers) main.cpp
	$(compiler) $(compilerFlags) -o main $(flags) $(objects) main.cpp $(paths) $(links)

# planar elements with WADG for heterogeneous media
main_subelem: $(objects) $(headers) main_subelem.cpp
	$(compiler) $(compilerFlags) -o main_subelem $(flags) $(objects) main_subelem.cpp $(paths) $(links)

# planar elements for piecewise constant media - can use Bernstein-Bezier
main_planar:	$(objects) $(headers) main_planar.cpp
	$(compiler) $(compilerFlags) -o main_planar $(flags) $(objects) main_planar.cpp $(paths) $(links)

$(objPath)/%.o:$(srcPath)/%.cpp $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(srcPath)/,$(incPath)/,$(<:.cpp=.tpp)))
	$(compiler) $(compilerFlags) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f $(objPath)/*.o;
	rm -f main;
	rm -f main_planar;
	rm -f main_subelem;
	rm -rf main.dSYM/;

.PHONY: main subelem planar clean
#=================================================
