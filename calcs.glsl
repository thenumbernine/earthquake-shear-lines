layout(local_size_x=<?=localSize.x
	?>, local_size_y=<?=localSize.y
	?>, local_size_z=<?=localSize.z
	?>) in;

<?=calcFlagsCode?>

real rad(real const d) {
	return d * M_PI / 180.;
}

constant const real WGS84_a = <?=glnumber(WGS84_a)?>.;		// m
constant const real WGS84_b = <?=glnumber(WGS84_b)?>;	// m @ polar radius
constant const real WGS84_esq = 1. - WGS84_b * WGS84_b / (WGS84_a * WGS84_a);
constant const real WGS84_e = <?=glnumber(WGS84_e)?>;//sqrt(WGS84_esq); ... but no constexpr in CL .c ... have to upgrade to clcpp but then you need to rope in spirv compiler ... too much drama ...
constant const real WGS84_flattening = 1. - WGS84_b / WGS84_a;
constant const real WGS84_inverseFlattening = 298.257223563;
constant const real WGS84_eccentricitySquared = (2. * WGS84_inverseFlattening - 1.) / (WGS84_inverseFlattening * WGS84_inverseFlattening);
constant const real gravitationalConstant = 6.6738480e-11;	// m^3 / (kg s^2)

real WGS84_calc_N(real const sinTheta) {
	real const denom = sqrt(1. - WGS84_eccentricitySquared * sinTheta * sinTheta);
	return WGS84_a / denom;
}

// how many copies of this do I need?
// This is the CL version to go with the GLSL version...
// I think I can gsub one to the other ...
real4 chart_WGS84(real4 const latLonHeight) {
	real const lat = latLonHeight.x;
	real const lon = latLonHeight.y;
	real const height = latLonHeight.z;

	real const phi = rad(lon);		// spherical φ
	real const theta = rad(lat);		// spherical inclination angle (not azumuthal θ)
	real const cosTheta = cos(theta);
	real const sinTheta = sin(theta);

	real const N = WGS84_calc_N(sinTheta);

	real const NPlusH = N + height;
	real4 y = (real4)(
		NPlusH * cosTheta * cos(phi),
		NPlusH * cosTheta * sin(phi),
		(N * (1. - WGS84_eccentricitySquared) + height) * sinTheta,
		0.
	);
	return y;
}

void chart_WGS84_basis(
	real4 const latLonHeight,
	real4 * const e_phi,
	real4 * const e_theta,
	real4 * const e_negR
) {
	real const phi = rad(latLonHeight.x);
	real const lambda = rad(latLonHeight.y);
	real const height = latLonHeight.z;

	real const cosLambda = cos(lambda);
	real const sinLambda = sin(lambda);

	real const cosPhi = cos(phi);
	real const sinPhi = sin(phi);
	real const dphi_cosPhi = -sinPhi;
	real const dphi_sinPhi = cosPhi;

	real const rCart = WGS84_a / sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	real const tmp = sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	real const dphi_rCart = WGS84_a / (tmp * tmp * tmp) * WGS84_esq * sinPhi * dphi_sinPhi;

	real const rCart_over_a = 1. / sqrt(1. - WGS84_esq * sinPhi * sinPhi);

	real const xp = (rCart + height) * cosPhi;
	real const dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi;
	real const dheight_xp = cosPhi;

	real const xp_over_a = (rCart_over_a + height / WGS84_a) * cosPhi;

	real const zp = (rCart * (1. - WGS84_esq) + height) * sinPhi;
	real const dphi_zp = (dphi_rCart * (1. - WGS84_esq)) * sinPhi + (rCart * (1. - WGS84_esq) + height) * dphi_sinPhi;
	real const dheight_zp = sinPhi;

	real const zp_over_a = (rCart_over_a * (1. - WGS84_esq) + height / WGS84_a) * sinPhi;

	real const r2D = sqrt(xp * xp + zp * zp);
	real const dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D;
	real const dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D;

	real const r2D_over_a = sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a);
	real const dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D;

	real const sinPhiSph = zp / r2D;
	real const dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D);
	real const dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D);

	real const cosPhiSph = sqrt(1. - sinPhiSph * sinPhiSph);
	//d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du;
	real const dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph;
	real const dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph;

	//real x = r2D * cosPhiSph / WGS84_a * cosLambda;
	//real y = r2D * cosPhiSph / WGS84_a * sinLambda;
	//real z = r2D * sinPhiSph / WGS84_a;

	real4 const dphi = (real4)(
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda,
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda,
		(dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph),
		0.);

	real4 const dlambda = (real4)(
		-sinLambda,
		cosLambda,
		0.,
		0.);

	real4 const dheight = (real4)(
		(dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * cosLambda,
		(dheight_r2D * cosPhiSph + r2D * dheight_cosPhiSph) * sinLambda,
		(dheight_r2D * sinPhiSph + r2D * dheight_sinPhiSph),
		0.);

	// are these normalized? dlambda looks like it
	*e_phi = dphi;		// d/d latlon.x
	*e_theta = dlambda;	// d/d latlon.y
	*e_negR = -dheight;	// d/d latlon.z s.t. it forms a rhs basis
}


//returns magnitude g+m-2d
real4 calcGravityAccel(
	real4 const posOnEarth, 			//relative to earth, in meters
	real4 const planetPosRelEarth 	//relative to earth, in meters, w = mass in kg
) {
	real4 x = posOnEarth - planetPosRelEarth;
	x.w = 0.;
	real const r = length(x);
	real const planetMass = planetPosRelEarth.w;
	return x * (-planetMass * gravitationalConstant / (r * r * r));
}

//returns magnitude g+m-2d
real4 calcTidalAccel(
	real4 const posOnEarth,
	real4 const normal,
	real4 const planetPosRelEarth
) {
	real4 x = posOnEarth - planetPosRelEarth;
	x.w = 0.;
	real const r = length(x);
	real const r2 = r * r;
	real const r3 = r * r2;
	real const xDotN = dot(x, normal);
	real const planetMass = planetPosRelEarth.w;
	return (x * (3. * xDotN / r2) - normal) * (gravitationalConstant * planetMass / r3);
}

real4 calcAccelAtPoint(
	real4 const posOnEarth,
	real4 const normal,
	global real4 const * const planetPosRelEarth,
	int const calcFlags
) {
	real4 accel = (real4)(0., 0., 0., 0.);

	for (int j = 0; j < <?=#Planets.planetClasses?>; ++j) {
	//for (int j = <?=Planets.indexes.earth-1?>; j <= <?=Planets.indexes.earth-1?>; ++j) {
		if (planetPosRelEarth[j].w > 0.) {
			if (calcFlags & calcFlags_calcTides) {
				real4 tideAccel = calcTidalAccel(posOnEarth, normal, planetPosRelEarth[j]);
				if ((calcFlags & calcFlags_calcTides) == calcFlags_calcTides) {
					accel += tideAccel;
				} else {
					real4 tideNormal = normal * dot(tideAccel, normal);
					real4 tideTangent = tideAccel - tideNormal;
					if (calcFlags & calcFlags_calcTidesVert) accel += tideNormal;
					if (calcFlags & calcFlags_calcTidesHorz) accel += tideTangent;
				}
			}
			if (calcFlags & calcFlags_calcGrav) {
				real4 gravAccel = calcGravityAccel(posOnEarth, planetPosRelEarth[j]);
				if ((calcFlags & calcFlags_calcGrav) == calcFlags_calcGrav) {
					accel += gravAccel;
				} else {
					real4 gravNormal = normal * dot(gravAccel, normal);
					real4 gravTangent = gravAccel - gravNormal;
					if (calcFlags & calcFlags_calcGravVert) accel += gravNormal;
					if (calcFlags & calcFlags_calcGravHorz) accel += gravTangent;
				}
			}
		}
	}

	if (calcFlags & calcFlags_calcCrossPos2D) {
		// why is this getting zeroes at the x and y axis ...
		accel = cross((real4)(posOnEarth.xy, 0., 0.), accel);	// torque = position cross accel
		// TODO what about velocity as well? going with vs against earth rotation?
	}
	if (calcFlags & calcFlags_calcCrossPos3D) {
		accel = cross(posOnEarth, accel);	// torque = position cross accel
		// TODO what about velocity as well? going with vs against earth rotation?
	}
	if (calcFlags & calcFlags_calcCrossVel) {
		real radsPerSec = 2. * M_PI / (1. - 1. / 365.25) / (60. * 60. * 24.);
		real4 velOnEarth = (real4)(
			radsPerSec * -posOnEarth.y,
			radsPerSec * posOnEarth.x,
			0.,
			0.
		);
		accel = cross(velOnEarth, accel);
	}

	return accel;
}

kernel void calcTideAndGrav(
	global real * const out,		// output[i] = force component specified
	global real4 const * const planetPosRelEarth,	// input positions[9] = (pos - earth.pos) x y z and mass = w
	int const calcFlags,
	int const display
) {
	initKernelForSize(<?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);
	// index, x = longitude = east/west, y = latitude = north/south

	// lat lon, x = latitude = north/south, longitude = left/right
	real4 latlon = (real4)(
		(real)((i.y + .5) / <?=glnumber(forceTexSize.y)?>) * 180. - 90.,
		(real)((i.x + .5) / <?=glnumber(forceTexSize.x)?>) * 360. - 180.,
		0.,
		0.
	);

	real4 const posOnEarth = chart_WGS84(latlon);

	real4 e_phi, e_theta, e_negR;
	chart_WGS84_basis(latlon, &e_phi, &e_theta, &e_negR);
	real4 const e_r = -normalize(e_negR);
	e_theta = normalize(e_theta);
	e_phi = normalize(e_phi);

	// should I be using e_r or normalize(posOnEarth) ? or are they the same?
	//real4 normal = normalize(posOnEarth);
	real4 const normal = e_r;

	//spherical coordinates, so theta goes from north pole to south pole, and phi goes around the earth
	real const rsurf = length(posOnEarth);
	real const sinTheta = length(posOnEarth.xy) / rsurf;
	real const dphi = <?=glnumber((2 * math.pi) / forceTexSize.x)?>;
	real const dtheta = <?=glnumber(math.pi / forceTexSize.y)?>;
	real const invdphi = <?=glnumber(forceTexSize.x / (2 * math.pi))?>;
	real const invdtheta = <?=glnumber(forceTexSize.y / math.pi)?>;

	real4 accel;
	if (calcFlags & calcFlags_calcRadialIntegral) {
		accel = (real4)(0., 0., 0., 0.);
		//if calcFlags_calcRadialIntegral then integrate ...
<?
local raddiv = 100
?>
		real const oneOverRadDiv = <?=glnumber(1 / raddiv)?>;
		real const dr = rsurf * oneOverRadDiv;
		real const sinTheta_dr_dtheta_dphi = sinTheta * dr * dtheta * dphi;
		for (int i = 0; i < <?=raddiv?>; ++i) {
			real const f = ((real)i + .5) * oneOverRadDiv;
			real const r = f * rsurf;
			real4 accelAtPt = calcAccelAtPoint(posOnEarth * f, normal, planetPosRelEarth, calcFlags);
			// TODO trapezoid integral across each face or something
			accel += accelAtPt * r * sinTheta_dr_dtheta_dphi;
		}
	} else {
		accel = calcAccelAtPoint(posOnEarth, normal, planetPosRelEarth, calcFlags);
	}

	// I'm breaking with chart input component convention and I'm putting r first here
	// so this is textbook math spherical basis convention of (r,theta,phi), not my geographic-chart convention of lat/lon/-height
	real4 const accel_r_theta_phi = (real4)(
		dot(accel, e_r),
		dot(accel, e_theta),
		dot(accel, e_phi),
		0.
	);

	// TODO display components, not always length
	// x y z r theta phi xy xz yz rth rph thph
	switch (display) {
// 1D:
	case displayComponent_x:
		out[index] = accel.x;
		break;
	case displayComponent_y:
		out[index] = accel.y;
		break;
	case displayComponent_z:
		out[index] = accel.z;
		break;
	case displayComponent_r:
		out[index] = accel_r_theta_phi.x;
		break;
	case displayComponent_theta:
		out[index] = accel_r_theta_phi.y;
		break;
	case displayComponent_phi:
		out[index] = accel_r_theta_phi.z;
		break;

// 2D:
	case displayComponent_angle_xy:
		out[index] = atan2(accel.y, accel.x);
		break;
	case displayComponent_angle_yz:
		out[index] = atan2(accel.z, accel.y);
		break;
	case displayComponent_angle_zx:
		out[index] = atan2(accel.x, accel.z);
		break;
	case displayComponent_length_xy:
		out[index] = length(accel.xy);
		break;
	case displayComponent_length_yz:
		out[index] = length(accel.yz);
		break;
	case displayComponent_length_zx:
		out[index] = length(accel.zx);
		break;

	case displayComponent_angle_rtheta:
		out[index] = atan2(accel_r_theta_phi.y, accel_r_theta_phi.x);
		break;
	case displayComponent_angle_thetaphi:
		out[index] = atan2(accel_r_theta_phi.z, accel_r_theta_phi.y);
		break;
	case displayComponent_angle_phir:
		out[index] = atan2(accel_r_theta_phi.x, accel_r_theta_phi.z);
		break;
	case displayComponent_length_rtheta:
		out[index] = length(accel_r_theta_phi.xy);
		break;
	case displayComponent_length_thetaphi:
		out[index] = length(accel_r_theta_phi.yz);
		break;
	case displayComponent_length_phir:
		out[index] = length(accel_r_theta_phi.zx);
		break;

// 3D:
	case displayComponent_length_xyz:
		out[index] = length(accel);
		break;
	default:
		out[index] = -(float)0xdeadbeef;
		break;
	}
}

//if (calcFlags & calcFlags_calcSurfaceGradMagn) {
kernel void gradient(
	global real * const out,
	global real const * const in
) {
	initKernelForSize(<?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);

	real4 const latlon = (real4)(
		(real)((i.y + .5) / <?=glnumber(forceTexSize.y)?>) * 180. - 90.,
		(real)((i.x + .5) / <?=glnumber(forceTexSize.x)?>) * 360. - 180.,
		0.,
		0.
	);

	real4 const posOnEarth = chart_WGS84(latlon);

	//spherical coordinates, so theta goes from north pole to south pole, and phi goes around the earth
	real const rsurf = length(posOnEarth);
	real const sinTheta = length(posOnEarth.xy) / rsurf;
	real const dphi = <?=glnumber((2 * math.pi) / forceTexSize.x)?>;
	real const dtheta = <?=glnumber(math.pi / forceTexSize.y)?>;
	real const invdphi = <?=glnumber(forceTexSize.x / (2 * math.pi))?>;
	real const invdtheta = <?=glnumber(forceTexSize.y / math.pi)?>;

	// spherical phi, not wgs84 phi ... smh
	int4 ixL = i;
	ixL.x = (ixL.x + <?=forceTexSize.x?> - 1) % <?=forceTexSize.x?>;
	int const index_xL = indexForInt4ForSize(ixL, <?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);

	int4 ixR = i;
	ixR.x = (ixR.x + 1) % <?=forceTexSize.x?>;
	int const index_xR = indexForInt4ForSize(ixR, <?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);

	int4 iyL = i;
	iyL.y--;
	if (iyL.y < 0) {
		iyL.y = 0;
		iyL.x = <?=forceTexSize.x?> - 1 - iyL.x;
	}
	int const index_yL = indexForInt4ForSize(iyL, <?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);

	int4 iyR = i;
	iyR.y++;
	if (iyR.y >= <?=forceTexSize.y?>) {
		iyR.y = <?=forceTexSize.y-1?>;
		iyR.x = <?=forceTexSize.x?> - 1 - iyR.x;
	}
	int const index_yR = indexForInt4ForSize(iyR, <?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);

	real const f_xL = in[index_xL];
	real const f_xR = in[index_xR];
	real const f_yL = in[index_yL];
	real const f_yR = in[index_yR];

	real const df_dtheta = -(f_yR - f_yL) * invdtheta;
	real const df_dphi = (f_xR - f_xL) * invdphi;
	real const rSq = rsurf * rsurf;
	real const g_theta_theta = rSq;
	real const g_phi_phi = rSq * sinTheta * sinTheta;
	out[index] = sqrt(
		df_dphi * df_dphi * g_theta_theta
		+ df_dtheta * df_dtheta * g_phi_phi
	);
}

kernel void rescale(
	write_only image2d_t outTex,	//complains invalid GL object (why?)
	global real const * const inv,
	real const minv,
	real const maxv
) {
	initKernelForSize(<?=forceTexSize.x?>, <?=forceTexSize.y?>, 1);
	float result = (float)((inv[index] - minv) / (maxv - minv));
	write_imagef(outTex, i.xy, (float4)(result, 0., 0., 0.));
}

// I woulda thought OpenCL would have this already
typedef char int8_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long int64_t;
typedef unsigned long uint64_t;

<?=smallBodyTypeCode?>

real4 rotateX(real4 const v, real const th) {
	real const c = cos(th);
	real const s = sin(th);
	return (real4)(
		v.x,
		v.y * c - v.z * s,
		v.y * s + v.z * c,
		0.);
}

real4 rotateZ(real4 const v, real const th) {
	real const c = cos(th);
	real const s = sin(th);
	return (real4)(
		v.x * c - v.y * s,
		v.x * s + v.y * c,
		v.z,
		0.);
}

// This is duplciated in basis.rua
realsb4 rotateFromSolarToEarthFrame(
	realsb4 v,
	real jday
) {
	return rotateZ(v,
		<?=glnumber(julianBaseAngleInRad)?>
		+ (fmod(jday, 1.) + .5) * -2. * M_PI * (
			// convert from sinodic day (24 hours = 360 degrees x ())
			// ... to sidereal day
			1. / (1. - 1. / 365.25)
			//1
			//(1 - 1 / 365.25)
		)
	);
}

//also in earthquakes.rua
realsb4 tiltFromSolarSystemToEarthFrame(realsb4 v) {
	realsb th = M_PI / 180. * ((23 + 1/60*(26 + 1/60*(21.4119))));
	return rotateX(v, th);
}

kernel void updateSmallBodies(
	global SmallBody * const bodies,
	realsb const julianDay,
	realsb4 const earthPosInSolarSystem
) {
	initKernelForSize(<?=numSmallBodies?>, 1, 1);

	global SmallBody * const ke = bodies + index;

	// wait does that mean resetDate should be whenever the A and B coeffs were generated?
	// TODO fix this in visualize-smallbodies too?
	// and TODO store it somewhere, maybe in the extra body_t_desc.lua file?
	realsb timeAdvanced = julianDay - <?=glnumber(smallbodies_julianResetDay)?>;

	int orbitType = ke->orbitType;

	//https://en.wikipedia.org/wiki/Kepler%27s_laws_of_planetary_motion#Position_as_a_function_of_time

	realsb meanAnomaly, meanMotion;
	if (orbitType == ORBIT_ELLIPTIC) {
		meanMotion = 2. * M_PI / ke->orbitalPeriod;
		meanAnomaly = ke->meanAnomalyAtEpoch + meanMotion * (julianDay - ke->epoch);
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		meanAnomaly = ke->meanAnomaly;
		meanMotion = ke->meanAnomaly / (julianDay - ke->timeOfPeriapsisCrossing);
	} else if (orbitType == ORBIT_PARABOLIC) {
		//error'got a parabolic orbit'
	} else {
		//error'here'
	}

	realsb eccentricity = ke->eccentricity;

	//solve eccentricAnomaly from meanAnomaly via Newton Rhapson
	//for elliptical orbits:
	//	f(E) = M - E + e sin E = 0
	//	f'(E) = -1 + e cos E
	//for parabolic oribts:
	//	f(E) = M - E - E^3 / 3
	//	f'(E) = -1 - E^2
	//for hyperbolic orbits:
	//	f(E) = M - e sinh(E) - E
	//	f'(E) = -1 - e cosh(E)
	realsb eccentricAnomaly = meanAnomaly;
	for (int i = 0; i < 10; ++i) {
		realsb func, deriv;
		if (orbitType == ORBIT_PARABOLIC) {
			func = meanAnomaly - eccentricAnomaly - eccentricAnomaly * eccentricAnomaly * eccentricAnomaly / 3.;
			deriv = -1. - eccentricAnomaly * eccentricAnomaly;
		} else if (orbitType == ORBIT_ELLIPTIC) {
			func = meanAnomaly - eccentricAnomaly + eccentricity * sin(eccentricAnomaly);
			deriv = -1. + eccentricity * cos(eccentricAnomaly);	//has zeroes ...
		} else if (orbitType == ORBIT_HYPERBOLIC) {
			func = meanAnomaly + eccentricAnomaly - eccentricity  * sinh(eccentricAnomaly);
			deriv = 1. - eccentricity * cosh(eccentricAnomaly);
		} else {
			//error'here'
		}

		realsb delta = func / deriv;
		if (fabs(delta) < 1e-15) break;
		eccentricAnomaly -= delta;
	}

	//TODO don't use meanMotion for hyperbolic orbits
	//realsb fractionOffset = timeAdvanced * meanMotion / (2 * M_PI);
	realsb theta = timeAdvanced * meanMotion;
	realsb pathEccentricAnomaly = eccentricAnomaly + theta;
	global realsb * A = ke->A;
	global realsb * B = ke->B;

	//matches above
	realsb dt_dE;
	realsb semiMajorAxisCubed = ke->semiMajorAxis * ke->semiMajorAxis * ke->semiMajorAxis;
	const realsb gravitationalParameter = <?=glnumber(smallbodies_gravitationalParameter)?>;
	if (orbitType == ORBIT_PARABOLIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (1. + pathEccentricAnomaly * pathEccentricAnomaly);
	} else if (orbitType == ORBIT_ELLIPTIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (1. - ke->eccentricity * cos(pathEccentricAnomaly));
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		dt_dE = sqrt(semiMajorAxisCubed / gravitationalParameter) * (ke->eccentricity * cosh(pathEccentricAnomaly) - 1.);
	}
	realsb dE_dt = 1. / dt_dE;
	realsb coeffA, coeffB;
	//realsb coeffDerivA, coeffDerivB;
	if (orbitType == ORBIT_PARABOLIC) {
		//...?
	} else if (orbitType == ORBIT_ELLIPTIC) {
		coeffA = cos(pathEccentricAnomaly) - ke->eccentricity;
		coeffB = sin(pathEccentricAnomaly);
		//coeffDerivA = -sin(pathEccentricAnomaly) * dE_dt;
		//coeffDerivB = cos(pathEccentricAnomaly) * dE_dt;
	} else if (orbitType == ORBIT_HYPERBOLIC) {
		coeffA = ke->eccentricity - cosh(pathEccentricAnomaly);
		coeffB = sinh(pathEccentricAnomaly);
		//coeffDerivA = -sinh(pathEccentricAnomaly) * dE_dt;
		//coeffDerivB = cosh(pathEccentricAnomaly) * dE_dt;
	}
	realsb posX = A[0] * coeffA + B[0] * coeffB;
	realsb posY = A[1] * coeffA + B[1] * coeffB;
	realsb posZ = A[2] * coeffA + B[2] * coeffB;
	//realsb velX = A[0] * coeffDerivA + B[0] * coeffDerivB;	//m/day
	//realsb velY = A[1] * coeffDerivA + B[1] * coeffDerivB;
	//realsb velZ = A[2] * coeffDerivA + B[2] * coeffDerivB;

	// TODO also rotate by julian day fraction?
	real4 pos = rotateFromSolarToEarthFrame(
		tiltFromSolarSystemToEarthFrame(
			(realsb4)(posX - earthPosInSolarSystem.x, posY - earthPosInSolarSystem.y, posZ - earthPosInSolarSystem.z, 0.)
		),
		julianDay
	);
	ke->pos[0] = pos.x;
	ke->pos[1] = pos.y;
	ke->pos[2] = pos.z;
	// I'm not using this ... but maybe I should be ...
	//ke->vel[0] = velX + parent.vel.s[0];
	//ke->vel[1] = velY + parent.vel.s[1];
	//ke->vel[2] = velZ + parent.vel.s[2];

//TODO CL/GL interop if you do use this
#if 0
	// now update buffers
	self.bodyToEarthArray[0+2*index].x = ke->pos[0];
	self.bodyToEarthArray[0+2*index].y = ke->pos[1];
	self.bodyToEarthArray[0+2*index].z = ke->pos[2];
	self.bodyToEarthArray[1+2*index].x = ke->pos[0];
	self.bodyToEarthArray[1+2*index].y = ke->pos[1];
	self.bodyToEarthArray[1+2*index].z = ke->pos[2];
#endif

	//ke == this.keplerianOrbitalElements;
	//ke->meanAnomaly = meanAnomaly;	// this doesn't change, does it?
	ke->eccentricAnomaly = eccentricAnomaly;	// this does change ... eccentricAnomaly and pos
	//ke->fractionOffset = fractionOffset;
	// but TODO if I'm always recalculating the pos & vel from the time0 = A & B parameter generation time (solarsystem/jpl-ssd-smallbody/parse.lua)
	//  then should I also recalculate from the initial eccentricAnomaly ?
	//  or is the initial eccentricAnomaly always zero?
}

