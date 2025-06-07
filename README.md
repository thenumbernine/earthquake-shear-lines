TO CALC:
- graph how many quakes have the sun or moon within 5 degrees ... tighten constraints and test all planets?
- calculate & graph how many quakes are on geodesics that also go through the sun or moon within .... 1 degree? 0.1 degree? idk?
- make a graph, per quake, of the falloff of "quakeAlignWithPreviousGeodesicAngleThreshold" - number of geodesics that align with this quake vs the threshold angle size.
	... because I've got it as tight as 0.1 degrees and every single quake still aligns with geodesics ... i.e. for every quake, there's 2 or more other quakes on Earth already forming a straight line with it.
- make a heatmap of *all geodesics* formed from 0.1 degrees & +/- 3 days

TO PROGRAM:
- use bisect method for esarch for time windows in the earthquakes and flares data
- when creating geodesics, filter out identical pairs to prevent those huge alpha stacking over high-quantity (low-magnitude) areas
	- make this optional
- click to select point and geodesic
- calibrate julianDay calculator
- should I just store all planet angles at all earthquake & flare times?
- S.H. basis for calcs?  If I add GBP data, that's already in SH coefficients ...
