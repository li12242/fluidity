FLFLAGS = -v2 -l
FLMLMODEL = square-convection
FLUIDITY = dfluidity
MESH = square-structured-linear
MODEL = square-convection

default: input

input: clean
	cp src/*.node src/*.ele src/*.edge .

clean: clean-mesh clean-run-debug 
	rm -fr *.ele *.node *.halo *.edge *checkpoint* square-convection.detectors.dat square-convection_* fluidity.[le]*
clean-mesh:
	rm -f $(MESH).edge $(MESH).ele $(MESH).node
clean-run:
	rm -f $(FLMLMODEL)_?*.vtu $(FLMLMODEL)_?*_checkpoint.edge $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.node $(FLMLMODEL)_?*_checkpoint.flml $(FLMLMODEL)_?*_?*_checkpoint.vtu $(FLMLMODEL).detectors $(FLMLMODEL).stat fluidity.err-0 fluidity.log-0 matrixdump matrixdump.info
	rm -rf *flredecomp*
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS) $(MODEL).flml
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS) $(MODEL).flml
