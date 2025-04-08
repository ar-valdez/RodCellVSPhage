#ifndef _CELLS_H_
#define _CELLS_H_

#include "tools.h"
#include "Array.h"
#include "Neighbours.h"
#include "Phages.h"

class UniformGrid;
struct Cell;

double getPressure(const Tensor T);

double getPressure(const DoubleCoord F);

void MoveCell(int cellID, UniformGrid& Grid, const Cell* old_cells, Cell* new_cells, \
              const int* NeighbourList, double dt, DoubleArray2D& Height, \
              CoordArray2D& Normal, DoubleArray2D& Wall, const int caseTest);

void GrowCell(Cell& cell, int cellID, double dt, int* dividingCells, int& numDivide, \
              EnvArray3D& Env, AgaArray2D** Wal, UniformGrid& Grid, bool Nutrients);

void DivideCell(int parentID, int daughterID, Cell* cells, UniformGrid& Grid, \
                const int* neighbours, DoubleArray2D& Wall, DoubleArray2D& Height, \
                CoordArray2D& Normal, double t);

bool CheckDivide(double DivCrit, double Lcurrent, double Loc_Ldiv);

#endif // _CELLS_H_
