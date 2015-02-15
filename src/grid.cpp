#include "grid.h"

Cell::~Cell() 
{
	delete[] nodesInd;
	delete[] edgesInd;
}

Edge::~Edge() 
{
	delete[] c;
}

Grid::~Grid() 
{
	delete[] nodes;
	delete[] cells;
	delete[] edges;
}

int Grid::findEdge(int n1, int n2)
{
	for (int iEdge = 0; iEdge < eCount; iEdge++)
	{
		if ((edges[iEdge].n1 == n1 && edges[iEdge].n2 == n2) || (edges[iEdge].n1 == n2 && edges[iEdge].n2 == n1))
		{
			return iEdge;
		}
	}
	return -1;
}

