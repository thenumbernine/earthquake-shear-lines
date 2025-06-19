#!/bin/sh
# make a movie centered on orbital plane, of 1 month, for all quakes, of all geodesics touching this quake (within 0.1 degrees) in a 3-day-prior window
./run.rua\
	width=1366\
	height=768\
	viewOrthoSize=1.4\
	weight_Equirectangular=0\
	weight_Azimuthal_equidistant=1\
	recenterOn=2\
	pickCirclesMethod=1\
	greatArcAlpha=.25\
	angleFade=10\
	drawFlares=false\
	forceTexAlpha=0\
	calcGravVert=false\
	lineWidth=2\
	latlonroll0convergeRate=1\
	makeMovie
