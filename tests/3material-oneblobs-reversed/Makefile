input: clean
	gmsh -2 src/channel.geo
	../../bin/gmsh2triangle --2d src/channel.msh

clean:
	rm -rf src/channel.msh *.stat *.vtu channel.node channel.ele channel.edge *.log-0 *.err-0 *checkpoint *.convergence *.png \
	matrixdump matrixdump.info
