#ifndef _PHAGES_H_
#define _PHAGES_H_

#include "tools.h"
#include "Array.h"
#include "Neighbours.h"

class UniformGrid;
struct Cell;

////////////////////////////////////////////////////////////////////////////////
// Math functions
////////////////////////////////////////////////////////////////////////////////

int RU_Numbers(int min, int max);

double RN_Numbers(double mu, double sigma);

double RLogN_Numbers(double mu, double sigma);

////////////////////////////////////////////////////////////////////////////////
// Phage specific functions
////////////////////////////////////////////////////////////////////////////////

void SourcePhage(int targetID, Cell* cells, double dt);

void BurstCell(int targetID, const int last_pos, Cell* cells, UniformGrid& Grid);

int UpdateInfection(int N_cells, Cell* cells, double dt);

int getAffectedList(int cellID, Cell* cells, UniformGrid& Grid, int maxNeighbours, int** NeighbourList, int* infectList);

void movePhage(double cx, double cy, double cz, double Phageupdate);

bool infectCell(DoubleCoord loc, int CellID, Cell* cells, double tolAds);

////////////////////////////////////////////////////////////////////////////////
// Properties of Cells
////////////////////////////////////////////////////////////////////////////////

double getCellRadius(int CellID, Cell* cells);

double getCellangle(int CellID, Cell* cells);

DoubleCoord getCellcenter(int CellID, Cell* cells);

int GetBoxNumber(UniformGrid& Grid, DoubleCoord loc);

void GetBoxCell(UniformGrid& Grid, DoubleCoord loc, int* CellinBox);

double getColonyRadii(int N_cells, Cell* cells);

double getColonyHeight(int N_cells, Cell* cells);

////////////////////////////////////////////////////////////////////////////////
// Here comes all the Cases PacMan, Border, Casper, etc
////////////////////////////////////////////////////////////////////////////////

int getBorderCell(int N_cells, Cell* cells, int srcs, int* CellPhageList);

int getCenterCell(int N_cells, Cell* cells, int srcs, int* CellPhageList);

int getPlaqueCell(int N_cells, Cell* cells, int srcs, int* PlaqueList);

int getPacmanCell(int N_cells, Cell* cells, int srcs, int* CellPhageList);




#endif // _PHAGES_H_
