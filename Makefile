lib.name = chaos

class.sources =  src/lorenz~.c  src/lorenzrk4~.c src/duffing~.c src/duffingeuler~.c src/doublependulum~.c src/rossler~.c src/lorenzattr~.c src/chua~.c src/chuaforce~.c src/henon~.c src/logistic~.c src/clogistic~.c src/tent~.c src/tentmap~.c src/rosslermod~.c

PDLIBBUILDER_DIR=pd-lib-builder/
include ${PDLIBBUILDER_DIR}/Makefile.pdlibbuilder
