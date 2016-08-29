PROJ_DIR:=$(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ifndef OCCA_DIR
  include $(PROJ_DIR)/../../../scripts/makefile
else
  include ${OCCA_DIR}/scripts/makefile
endif

#flags += -DtFloat=float -DOCCA_GL_ENABLED=1 -Dp_N=$(N) -DNDG3d -g
paths += -I./include

B ?= 0 # if B not set, default to B = 0 (nodal basis)
flags += -DOCCA_GL_ENABLED=1 -Dp_N=$(N) -DUSE_BERN=$(B) -g
#flags += -DOCCA_GL_ENABLED=1 -Dp_N=$(N) -g

ifeq ($(OS),OSX)
	links += -framework OpenGL -framework GLUT
endif

ifeq ($(OS),LINUX)
	links +=
# -lGLU -lglut
endif


#---[ COMPILATION ]-------------------------------
headers = $(wildcard $(iPath)/*.hpp) $(wildcard $(iPath)/*.tpp)
sources = $(wildcard $(sPath)/*.cpp)

objects  = $(subst $(sPath)/,$(oPath)/,$(sources:.cpp=.o))

# curved elements with WADG (can also add heterogeneous)
main: $(objects) $(headers) main.cpp
	$(compiler) $(compilerFlags) -o main $(flags) $(objects) main.cpp $(paths) $(links)

# planar elements with WADG for heterogeneous media
subelem: $(objects) $(headers) main_subelem.cpp
	$(compiler) $(compilerFlags) -o main_subelem $(flags) $(objects) main_subelem.cpp $(paths) $(links)

# planar elements for piecewise constant media - can use Bernstein-Bezier
planar:	$(objects) $(headers) main_planar.cpp
	$(compiler) $(compilerFlags) -o main_planar $(flags) $(objects) main_planar.cpp $(paths) $(links)

$(oPath)/%.o:$(sPath)/%.cpp $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.hpp))) $(wildcard $(subst $(sPath)/,$(iPath)/,$(<:.cpp=.tpp)))
	$(compiler) $(compilerFlags) -o $@ $(flags) -c $(paths) $<

clean:
	rm -f $(oPath)/*.o;
	rm -f main;
	rm -f main_planar;
	rm -f main_subelem;
	rm -rf main.dSYM/;

.PHONY: main main_subelem clean
#=================================================
