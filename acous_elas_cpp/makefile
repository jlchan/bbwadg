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
headers = $(wildcard $(incPath)/*.hpp) $(wildcard $(incPath)/*.tpp)
sources = $(wildcard $(srcPath)/*.cpp)

objects  = $(subst $(srcPath)/,$(objPath)/,$(sources:.cpp=.o))

# curved elements with WADG (can also add heterogeneous)
main: $(objects) $(headers) main.cpp
	$(compiler) $(compilerFlags) -o main $(flags) $(objects) main.cpp $(paths) $(links)

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
