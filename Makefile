all:
	gcc rpoly.c MatLabUtilites.c MatLabUtilitesTest.c -o   matlabutilites.out -lm
	./matlabutilites.out
	rm matlabutilites.out
	gcc rpoly.c MatLabUtilites.c unit.c doubler.c angl.c newtonnu.c PreliminaryOrbitDeterminationTest.c -o   PreliminaryOrbitDeterminationTest.out -lm
	./PreliminaryOrbitDeterminationTest.out
	rm PreliminaryOrbitDeterminationTest.out