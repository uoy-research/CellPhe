// SaveInterpolatedRois.ijm is an ImageJ macro designed for interpolation
// of ROIs output by TrackMate.

// To run the macro, ensure that the ROI set output from TrackMate is loaded
// in the FIJI ROI manager. Change the filepath in the macro below to save ROIs 
// to your choice of directory prior to running.

count = roiManager("count");

for (i = 0; i < count; i++) {
	roiManager("Select", i);
	run("Interpolate", "interval=1");
	roiManager("update");
}
roiManager("deselect");
filepath = "PATH/FOR/OUTPUT/ROIS.ZIP"
roiManager("save", filepath);
