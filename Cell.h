/*
 * Cell.h
 *
 *  Created on: Dec 15, 2011
 *      Author: mya
 */

#ifndef CELL_H_
#define CELL_H_

#include "tools.h"
#include "Constants.h"


// Cell structure contains information specific to each cell
struct Cell
{
	Segment Position; 	// a segment has two coordinates p and q for the two vertices of the spherocylinder
	double Length;
	double Radius;
	double GrowthRate = maxGrowthRate;
	DoubleCoord Velocity;
	DoubleCoord AngularVelocity; // now a vector
	int Type;
	int Ancestor;
	int Nb           = 1;
	double PhageCell = 0.0;
	double AgeCell   = 0.0;
	double Ldiv      = L_divide;
	double Strain    = 0.0;
	double Pressure  = 0.0;
	bool Abortive    = false;
	bool Shrink      = false;
	bool Short       = false;
	bool Resistant   = false;
};

#endif /* CELL_H_ */
