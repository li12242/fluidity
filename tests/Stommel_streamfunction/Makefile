input:
	gmsh -2 -o Stommel_gyre.msh src/Stommel_gyre.geo
	../../bin/gmsh2triangle --2d Stommel_gyre.msh
	gmsh -2 -o Stommel_gyre_A.msh src/Stommel_gyre_A.geo
	../../bin/gmsh2triangle --2d Stommel_gyre_A.msh


clean:
	rm -f *.node *.ele *.edge *.msh *.vtu *.log* *.err* *.stat *dump* \
	matrixdump matrixdump.info
