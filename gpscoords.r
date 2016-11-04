#!/usr/bin/R	

# GPS coordinates-related functions


geodetic.distance <- function(point1, point2) { 
	# taken from http://www.biostat.umn.edu/~sudiptob/Software/distonearth.R
	R = 6371 
	p1rad = point1 * pi/180 
	p2rad = point2 * pi/180 
	d = sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))	
	d = acos(d) 
	return(R*d) 
} 

minutesseconds2decimal.coords = function(dms){
	degrees = dms[1]
	minutes = dms[2]
	seconds = dms[3]
	return(degrees + minutes/60 + seconds/3600)
}

gps2xyz.coords = function(lon, lat=NULL, radians=F){
	# takes in 2 numeric args, or a 2-vector as a single arg, describing longitude and latitude
	if (is.null(lat)){
		gps.coord = lon
		lon = gps.coord[1]
		lat = gps.coord[2]
	}
	if (!radians){
		lon = lon / 180
		lat = lat / 180
	}
	x = cos(lat) * cos(lon)
	y = cos(lat) * sin(lon)
	z = sin(lat)
	return(c(x,y,z))
}

xyz2gps.coords = function(x, y=NULL, z=NULL, radians=F){
	if (is.null(y) | is.null(z)){
		xyz.coord = x
		x = xyz.coord[1]
		y = xyz.coord[2]
		z = xyz.coord[3]
	}
	lon = atan2(y, x)
	hyp = sqrt(x^2 + y^2)
	lat = atan2(z, hyp)
	if (!radians){
		lon = lon * 180
		lat = lat * 180
	}
	return(c(lon, lat))
}

average.gps.coords = function(gps.coords, ...){
	# takes in a list of 2-vectors or data.frame describing each point's longitude and latitude
	if (is.data.frame(gps.coords)){
		xyz.coords = t(apply(gps.coords, 1, gps2xyz.coords, ...))
	}else{ if (is.list(gps.coords)){
		xyz.coords = t(sapply(gps.coords, gps2xyz.coords, ...))
	}}
	ave.xyz.coord = apply(xyz.coords, 2, mean)
	ave.gps.coord = xyz2gps.coords(ave.xyz.coord, ...)
	return(ave.gps.coord)
}
