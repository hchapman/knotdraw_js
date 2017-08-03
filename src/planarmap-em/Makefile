CC = emcc
OPT =

OBJS = PMdef.bc PMconjugation.bc PMenlight.bc PMextract.bc PMdisplay.bc \
       PMplanmap.bc interface.bc stats.bc PMemscripten.bc

all: $(OBJS) planarmap

PMdef.bc: PMdef.c
	$(CC) $(OPT) -c PMdef.c -o PMdef.bc

PMconjugation.bc: PMconjugation.c
	$(CC) $(OPT) -c PMconjugation.c -o PMconjugation.bc

PMenlight.bc: PMenlight.c PMenlight.h
	$(CC) $(OPT) -c PMenlight.c -o PMenlight.bc

PMextract.bc: PMextract.c
	$(CC) $(OPT) -c PMextract.c -o PMextract.bc

PMdisplay.bc: PMdisplay.c PMdisplay.h
	$(CC) $(OPT) -c PMdisplay.c -o PMdisplay.bc

PMplanmap.bc: PMplanmap.c PMplanmap.h
	$(CC) $(OPT) -c PMplanmap.c -o PMplanmap.bc

PMemscripten.bc: PMemscripten.c PMemscripten.h
	$(CC) $(OPT) -c PMemscripten.c -o PMemscripten.bc

stats.bc: stats.c stats.h
	$(CC) $(OPT) -c stats.c -o stats.bc

planarmap: $(OBJS)
	$(CC) $(OPT) $(OBJS) -lm -o libplanarmap-em.js

interface.bc: interface.c
	$(CC) $(OPT) -c interface.c -o interface.bc

clean:
	/bin/rm -f $(OBJS)

VERSION=PlanarMap-v1-2b
