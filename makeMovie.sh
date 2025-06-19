#!/bin/sh
# make a movie of 1 month, for all quakes, of all geodesics touching any quake (within 0.1 degrees) in a 3-day-prior window
./run.rua\
	width=768\
	height=768\
	viewOrthoSize=1.4\
	weight_Equirectangular=0\
	weight_Azimuthal_equidistant=1\
	pickCirclesMethod=2\
	greatArcAlpha=.25\
	angleFade=10\
	drawFlares=false\
	forceTexAlpha=0\
	calcGravVert=false\
	makeMovie
