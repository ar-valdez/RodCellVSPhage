#ifndef INTEGRATE_H_
#define INTEGRATE_H_
#include "UniformGrid.h"
#include "tools.h"
#include "Array.h"
#include "Cell.h"

void integrate(const double dt, int cellID, const Cell* old_cells, Cell* new_cells, \
               const int* neighbours, DoubleArray2D& Height, CoordArray2D& Normal, \
               UniformGrid& Grid, const IntCoord& XYAddress, DoubleArray2D& Wall,
               const int caseTest);

void getPositions(const double dt, const DoubleCoord& F, const DoubleCoord& T, const Cell& old_cell, Cell& new_cell);

#endif /* INTEGRATE_H_ */
