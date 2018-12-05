all: bench

LIKWID_INCDIR=/apps/likwid/4.1.1/include
LIKWID_LIBDIR=apps/likwid/4.1.1/lib

bench: bench.c
	gcc -O3 -fopenmp -DLIKWID_PERFMON -I$LIKWID_INCDIR -L$LIKWID_LIBDIR bench.c -o stream_likwid.exe -llikwid
