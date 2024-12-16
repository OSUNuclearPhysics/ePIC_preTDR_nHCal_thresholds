/*
 * BasicUtil.h
 *
 *  Created on: 14 gru 2023
 *      Author: Khaless
 */

#ifndef BASICUTIL_H_
#define BASICUTIL_H_


#include <TMath.h>
#include "TVector.h"


double *calculateDirection(double eta, double phi, double r = 1.0);
TVector3 calculateDirectionVec(double eta, double phi, double r = 1.0);

TVector3 projTrack(double eta, double phi, double r = 80000.0, TVector3 start = TVector3(0.0, 0.0, 0.0)); // eta, phi, rmax[mm]
TVector3 projTrackZ(double eta, double phi, double z = -3950.0, TVector3 start = TVector3(0.0, 0.0, 0.0)); // eta, phi, zmax[mm]



double *calculateDirection(double eta, double phi, double r)
{
	double *dir = new double[3];

	double theta = 2.0*atan(exp(-eta));
	//double theta = 2.0*TMath::ATan(TMath::Exp(-eta));

	TVector3 vec;
	vec.SetMag(r);
	vec.SetTheta(theta);
	vec.SetPhi(phi);

	dir[0] = vec.x();
	dir[1] = vec.y();
	dir[2] = vec.z();

	return dir;
}


TVector3 calculateDirectionVec(double eta, double phi, double r)
{
	double theta = 2.0*atan(exp(-eta));
	//double theta = 2.0*TMath::ATan(TMath::Exp(-eta));

	TVector3 vec(1,1,1);
	vec.SetMag(r);
	vec.SetTheta(theta);
	vec.SetPhi(phi);

	return vec;
}


TVector3 projTrack(double eta, double phi, double r, TVector3 start)
{

	 TVector3 dir = calculateDirectionVec(eta, phi);

	 TVector3 endPointVec(start.x(), start.y(), start.z());
	 TVector3 diffPointVec(0.0, 0.0, 0.0);
/*
	 diffPointVec.SetMag(r);
	 diffPointVec.SetTheta(dir.Theta());
	 diffPointVec.SetPhi(phi);

	 diffPointVec=diffPointVec-start;
*/

	 endPointVec = start+diffPointVec;

	 return endPointVec;

}

TVector3 projTrackZ(double eta, double phi, double z, TVector3 start)
{

	 TVector3 dir = calculateDirectionVec(eta, phi);

	 TVector3 endPointVec(start.x(), start.y(), start.z());
	 TVector3 diffPointVec(0.0, 0.0, z-start.z());

	 double deltaZ = diffPointVec.z();
	 double deltaXY = deltaZ/(sinh(eta));

	 diffPointVec.SetXYZ(deltaXY*cos(phi), deltaXY*sin(phi), deltaZ);

	 endPointVec = start+diffPointVec;

	 return endPointVec;

}

#endif /* BASICUTIL_H_ */
