count = roiManager("count");

for (i = 0; i < count; i++) {
	roiManager("Select", i);
	run("Interpolate", "interval=1");
	roiManager("update");
}
Experiment = "21062022_C3_4_Red";
roiManager("deselect");
filepath = "/Users/laurawiggins/Desktop/" + Experiment + "-newROIs.zip";
roiManager("save", filepath);