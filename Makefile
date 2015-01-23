eBpFun.o: eBpFun.c head/constant.h head/eBpFun.h head/rhoFun.h
	gcc -c eBpFun.c -O3

eBsFun.o: eBsFun.c head/constant.h head/eBsFun.h head/rhoFun.h
	gcc -c eBsFun.c -O3

rhoFun.o: rhoFun.c head/rhoFun.h
	gcc -c rhoFun.c -O3

eB.o: eB.c head/eB.h head/eBpFun.h head/eBsFun.h
	gcc -c eB.c -O3

CME.o: CME.c head/eB.h
	gcc -c CME.c -O3

findMax.o: findMax.c head/eB.h
	gcc -c findMax.c -O3

CME: CME.o eB.o rhoFun.o eBpFun.o eBsFun.o 
	gcc CME.o eB.o rhoFun.o eBpFun.o eBsFun.o -o CME -lgsl -lgslcblas -lm -O3


findMax: findMax.o eB.o rhoFun.o eBpFun.o eBsFun.o
	gcc findMax.o eB.o rhoFun.o eBpFun.o eBsFun.o -o findMax -lgsl -lgslcblas -lm -O3
