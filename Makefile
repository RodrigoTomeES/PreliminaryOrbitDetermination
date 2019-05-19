all:
	gcc rpoly.c MatLabUtilites.c MatLabUtilitesTest.c -o   matlabutilites.out -lm
	./matlabutilites.out

	gcc rpoly.c MatLabUtilites.c unit.c doubler.c angl.c newtonnu.c R_x.c R_y.c R_z.c Frac.c timediff.c Mjday.c Position.c MeanObliquity.c NutAngles.c gibbs.c gmst.c EqnEquinox.c IERS.c gast.c NutMatrix.c PoleMatrix.c PrecMatrix.c GHAMatrix.c lambert_gooding.c rv2coe.c hgibbs.c anglesdr.c anglesg.c PreliminaryOrbitDeterminationTest.c -o PreliminaryOrbitDeterminationTest.out -lm
	./PreliminaryOrbitDeterminationTest.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesdr.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gmst.c EqnEquinox.c Frac.c example1.c -o example1.out -lm
	./example1.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesdr.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gmst.c EqnEquinox.c Frac.c example2.c -o example2.out -lm
	./example2.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesdr.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gmst.c EqnEquinox.c Frac.c example3.c -o example3.out -lm
	./example3.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesg.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gibbs.c hgibbs.c rv2coe.c angl.c newtonnu.c example5.c -o example5.out gmst.c EqnEquinox.c Frac.c -lm
	./example5.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesdr.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gmst.c EqnEquinox.c Frac.c example6.c -o example6.out -lm
	./example6.out

	gcc MatLabUtilites.c Mjday.c Position.c rpoly.c timediff.c IERS.c PrecMatrix.c NutMatrix.c PoleMatrix.c GHAMatrix.c anglesdr.c R_z.c R_y.c lambert_gooding.c MeanObliquity.c NutAngles.c R_x.c gast.c doubler.c unit.c gmst.c EqnEquinox.c Frac.c example7.c -o example7.out -lm
	./example7.out
