FLFLAGS = -v2 -l
FLMLMODEL = stokes-1e4-p2p1-parallel
FLUIDITY = fluidity
MESH = square-structured-linear
MODEL = Stokes-square-convection-1e4-parallel-p2p1

default: input

input: src/$(MESH).edge src/$(MESH).ele src/$(MESH).node clean
	cp src/$(MESH).edge src/$(MESH).ele src/$(MESH).node .

clean: clean-mesh clean-run-debug
clean-mesh:
	rm -f $(MESH).ele $(MESH).face $(MESH).node $(MESH).edge
	rm -f $(MESH)_?*.ele $(MESH)_?*.edge $(MESH)_?*.face $(MESH)_?*.halo $(MESH)_?*.node
clean-run:
	rm -f $(FLMLMODEL)_?*.pvtu $(FLMLMODEL)_?*.vtu
	rm -f $(FLMLMODEL).detectors $(FLMLMODEL).stat tmpf*
	rm -f fluidity.err-?* fluidity.log-?*
	rm -f matrixdump matrixdump.info 
	rm -f $(FLMLMODEL)_?*_checkpoint.flml
	rm -f $(FLMLMODEL)_?*_checkpoint_?*.ele $(FLMLMODEL)_?*_checkpoint_?*.face $(FLMLMODEL)_?*_checkpoint_?*.halo $(FLMLMODEL)_?*_checkpoint_?*.node
	rm -f $(FLMLMODEL)_?*_?*_checkpoint.pvtu $(FLMLMODEL)_?*_?*_checkpoint.vtu
	rm -f $(FLMLMODEL)_?*_checkpoint.ele $(FLMLMODEL)_?*_checkpoint.face $(FLMLMODEL)_?*_checkpoint.node
	rm -f $(FLMLMODEL)_?*_checkpoint.vtu
	rm -f -R $(FLMLMODEL)_* rm $(FLMLMODEL)*.dat
	rm -rf *flredecomp*
clean-run-debug: clean-run
	rm -f adapted_mesh_?*.vtu bounding_box_?*_?*.vtu final_metric_?*.vtu gmon.out gradation_metric_?*.vtu interpolation_metric_final_?*.vtu interpolation_metric_hessian_?*_?*.vtu interpolation_metric_merge_?*_?*.vtu interpolation_metric_metric_?*_?*.vtu metric_input_?*_?*.vtu

run: clean-run input
	$(FLUIDITY) $(FLFLAGS) $(MODEL).flml $MESH.*
run-debug: clean-run-debug input
	$(FLUIDITY)-debug $(FLFLAGS) $(MODEL).flml


