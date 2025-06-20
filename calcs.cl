<?=calcFlagsCode?>

real rad(real d) {
	return d * M_PI / 180.;
}

constant const real WGS84_a = 6378137.;		// m
constant const real WGS84_b = 6356752.3142;	// m @ polar radius
constant const real WGS84_esq = 1. - WGS84_b * WGS84_b / (WGS84_a * WGS84_a);
constant const real WGS84_e = <?=clnumber(WGS84_e)?>;//sqrt(WGS84_esq); ... but no constexpr in CL .c ... have to upgrade to clcpp but then you need to rope in spirv compiler ... too much drama ...
constant const real WGS84_flattening = 1. - WGS84_b / WGS84_a;
constant const real WGS84_inverseFlattening = 298.257223563;
constant const real WGS84_eccentricitySquared = (2. * WGS84_inverseFlattening - 1.) / (WGS84_inverseFlattening * WGS84_inverseFlattening);
constant const real gravitationalConstant = 6.6738480e-11;	// m^3 / (kg s^2)

real WGS84_calc_N(real sinTheta) {
	real denom = sqrt(1. - WGS84_eccentricitySquared * sinTheta * sinTheta);
	return WGS84_a / denom;
}

// how many copies of this do I need?
// This is the CL version to go with the GLSL version...
// I think I can gsub one to the other ...
real4 chart_WGS84(real4 latLonHeight) {
	real lat = latLonHeight.x;
	real lon = latLonHeight.y;
	real height = latLonHeight.z;

	real phi = rad(lon);		// spherical φ
	real theta = rad(lat);		// spherical inclination angle (not azumuthal θ)
	real cosTheta = cos(theta);
	real sinTheta = sin(theta);

	real N = WGS84_calc_N(sinTheta);

	real NPlusH = N + height;
	real4 y = (real4)(
		NPlusH * cosTheta * cos(phi),
		NPlusH * cosTheta * sin(phi),
		(N * (1. - WGS84_eccentricitySquared) + height) * sinTheta,
		0.
	);
	return y;
}

void chart_WGS84_basis(
	real4 latLonHeight,
	real4 *e_phi,
	real4 *e_theta,
	real4 *e_negR
) {
	real phi = rad(latLonHeight.x);
	real lambda = rad(latLonHeight.y);
	real height = latLonHeight.z;

	real cosLambda = cos(lambda);
	real sinLambda = sin(lambda);

	real cosPhi = cos(phi);
	real sinPhi = sin(phi);
	real dphi_cosPhi = -sinPhi;
	real dphi_sinPhi = cosPhi;

	real rCart = WGS84_a / sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	real tmp = sqrt(1. - WGS84_esq * sinPhi * sinPhi);
	real dphi_rCart = WGS84_a / (tmp * tmp * tmp) * WGS84_esq * sinPhi * dphi_sinPhi;

	real rCart_over_a = 1. / sqrt(1. - WGS84_esq * sinPhi * sinPhi);

	real xp = (rCart + height) * cosPhi;
	real dphi_xp = dphi_rCart * cosPhi + (rCart + height) * dphi_cosPhi;
	real dheight_xp = cosPhi;

	real xp_over_a = (rCart_over_a + height / WGS84_a) * cosPhi;

	real zp = (rCart * (1. - WGS84_esq) + height) * sinPhi;
	real dphi_zp = (dphi_rCart * (1. - WGS84_esq)) * sinPhi + (rCart * (1. - WGS84_esq) + height) * dphi_sinPhi;
	real dheight_zp = sinPhi;

	real zp_over_a = (rCart_over_a * (1. - WGS84_esq) + height / WGS84_a) * sinPhi;

	real r2D = sqrt(xp * xp + zp * zp);
	real dphi_r2D = (xp * dphi_xp + zp * dphi_zp) / r2D;
	real dheight_r2D = (xp * dheight_xp + zp * dheight_zp) / r2D;

	real r2D_over_a = sqrt(xp_over_a * xp_over_a + zp_over_a * zp_over_a);
	real dphi_r2D_over_a = (xp_over_a * dphi_xp + zp_over_a * dphi_zp) / r2D;

	real sinPhiSph = zp / r2D;
	real dphi_sinPhiSph = (dphi_zp * r2D - zp * dphi_r2D) / (r2D * r2D);
	real dheight_sinPhiSph = (dheight_zp * r2D - zp * dheight_r2D) / (r2D * r2D);

	real cosPhiSph = sqrt(1. - sinPhiSph * sinPhiSph);
	//d/du sqrt(1 - x^2) = -x/sqrt(1 - x^2) dx/du;
	real dphi_cosPhiSph = -sinPhi / cosPhiSph * dphi_sinPhiSph;
	real dheight_cosPhiSph = -sinPhi / cosPhiSph * dheight_sinPhiSph;

	//real x = r2D * cosPhiSph / WGS84_a * cosLambda;
	//real y = r2D * cosPhiSph / WGS84_a * sinLambda;
	//real z = r2D * sinPhiSph / WGS84_a;

	real4 dphi = (real4)(
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * cosLambda,
		(dphi_r2D_over_a * cosPhiSph + r2D_over_a * dphi_cosPhiSph) * sinLambda,
		(dphi_r2D_over_a * sinPhiSph + r2D_over_a * dphi_sinPhiSph),
		0.);

	real4 dlambda = (real4)(
		-sinLambda,
		cosLambda,
		0.,
		0.);

	real4 dheight = (real4)(
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
	real4 posOnEarth, 			//relative to earth, in meters
	real4 planetPosRelEarth 	//relative to earth, in meters, w = mass in kg
) {
	real4 x = posOnEarth - planetPosRelEarth;
	x.w = 0.;
	real const r = length(x);
	real const planetMass = planetPosRelEarth.w;
	return x * (-planetMass * gravitationalConstant / (r * r * r));
}

//returns magnitude g+m-2d
real4 calcTidalAccel(
	real4 posOnEarth,
	real4 normal,
	real4 planetPosRelEarth
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
	real4 posOnEarth,
	real4 normal,
	global real4 * const planetPosRelEarth,
	int calcFlags
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
	initKernelForSize(<?=clsize[1]?>, <?=clsize[2]?>, 1);
	// index, x = longitude = east/west, y = latitude = north/south

	// lat lon, x = latitude = north/south, longitude = left/right
	real4 latlon = (real4)(
		(real)((i.y + .5) / <?=clnumber(clsize[2])?>) * 180. - 90.,
		(real)((i.x + .5) / <?=clnumber(clsize[1])?>) * 360. - 180.,
		0.,
		0.
	);

	real4 posOnEarth = chart_WGS84(latlon);

	real4 e_phi, e_theta, e_negR;
	chart_WGS84_basis(latlon, &e_phi, &e_theta, &e_negR);
	real4 e_r = -normalize(e_negR);
	e_theta = normalize(e_theta);
	e_phi = normalize(e_phi);

	// should I be using e_r or normalize(posOnEarth) ? or are they the same?
	//real4 normal = normalize(posOnEarth);
	real4 normal = e_r;

	//spherical coordinates, so theta goes from north pole to south pole, and phi goes around the earth
	real const rsurf = length(posOnEarth);
	real const sinTheta = length(posOnEarth.xy) / rsurf;
	real const dphi = <?=clnumber((2 * math.pi) / clsize[1])?>;
	real const dtheta = <?=clnumber(math.pi / clsize[2])?>;
	real const invdphi = <?=clnumber(clsize[1] / (2 * math.pi))?>;
	real const invdtheta = <?=clnumber(clsize[2] / math.pi)?>;

	real4 accel;
	if (calcFlags & calcFlags_calcRadialIntegral) {
		accel = (real4)(0., 0., 0., 0.);
		//if calcFlags_calcRadialIntegral then integrate ...
<?
local raddiv = 100
?>
		real const oneOverRadDiv = <?=clnumber(1 / raddiv)?>;
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
	real4 accel_r_theta_phi = (real4)(
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
	initKernelForSize(<?=clsize[1]?>, <?=clsize[2]?>, 1);

	real4 latlon = (real4)(
		(real)((i.y + .5) / <?=clnumber(clsize[2])?>) * 180. - 90.,
		(real)((i.x + .5) / <?=clnumber(clsize[1])?>) * 360. - 180.,
		0.,
		0.
	);

	real4 posOnEarth = chart_WGS84(latlon);

	//spherical coordinates, so theta goes from north pole to south pole, and phi goes around the earth
	real const rsurf = length(posOnEarth);
	real const sinTheta = length(posOnEarth.xy) / rsurf;
	real const dphi = <?=clnumber((2 * math.pi) / clsize[1])?>;
	real const dtheta = <?=clnumber(math.pi / clsize[2])?>;
	real const invdphi = <?=clnumber(clsize[1] / (2 * math.pi))?>;
	real const invdtheta = <?=clnumber(clsize[2] / math.pi)?>;

	// spherical phi, not wgs84 phi ... smh
	int4 ixL = i;
	ixL.x = (ixL.x + <?=clsize[1]?> - 1) % <?=clsize[1]?>;
	int const index_xL = indexForInt4ForSize(ixL, <?=clsize[1]?>, <?=clsize[2]?>, 1);
	
	int4 ixR = i;
	ixR.x = (ixR.x + 1) % <?=clsize[1]?>;
	int const index_xR = indexForInt4ForSize(ixR, <?=clsize[1]?>, <?=clsize[2]?>, 1);
	
	int4 iyL = i;
	iyL.y--;
	if (iyL.y < 0) {
		iyL.y = 0;
		iyL.x = <?=clsize[1]?> - 1 - iyL.x;
	}
	int const index_yL = indexForInt4ForSize(iyL, <?=clsize[1]?>, <?=clsize[2]?>, 1);
	
	int4 iyR = i;
	iyR.y++;
	if (iyR.y >= <?=clsize[2]?>) {
		iyR.y = <?=clsize[2]-1?>;
		iyR.x = <?=clsize[1]?> - 1 - iyR.x;
	}
	int const index_yR = indexForInt4ForSize(iyR, <?=clsize[1]?>, <?=clsize[2]?>, 1);

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
<? if useGLSharing then ?>
	write_only image2d_t outTex,
<? else ?>
	global float * const outv,
<? end ?>
	global real const * const inv,
	real const minv,
	real const maxv
) {
	initKernelForSize(<?=clsize[1]?>, <?=clsize[2]?>, 1);
	float result = (float)((inv[index] - minv) / (maxv - minv));
<? if useGLSharing then ?>
	write_imagef(outTex, i.xy, (float4)(result, 0., 0., 0.));
<? else ?>
	outv[index] = result;
<? end ?>
}
