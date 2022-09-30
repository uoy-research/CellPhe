count = roiManager("count");

for (i = 0; i < count; i++) {
	roiManager("Select", i);
	run("Interpolate", "interval=1");
	roiManager("update");
}
roiManager("deselect");
filepath = "PATH/FOR/OUTPUT/ROIS.ZIP"
roiManager("save", filepath);
