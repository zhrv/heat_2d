#include "fvm_heat.h"
#include "tinyxml.h"
#include <string>
#include "global.h"

void FVM_Heat::init(char * xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		log("ERROR: %s\n", doc.ErrorDesc());
		exit(doc.ErrorId());
	}
	
	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );


	node0 = task->FirstChild("control");
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("CFL")->ToElement()->Attribute("value", &CFL);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	// чтение параметров о МАТЕРИАЛАХ
	node0 = task->FirstChild("materials");
	TiXmlNode* matNode = node0->FirstChild("material"); 
	while (matNode != NULL)
	{
		Material mat;
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		materials.push_back(mat);

		matNode = matNode->NextSibling("material");
	}
	matCount = materials.size();

	// чтение параметров о РЕГИОНАХ
	node0 = task->FirstChild("regions");
	TiXmlNode* regNode = node0->FirstChild("region");
	while (regNode != NULL)
	{
		Region reg;
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);
		
		node1 = regNode->FirstChild("parameters");
		node1->FirstChild( "T"  )->ToElement()->Attribute( "value", &reg.par.T );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &reg.par.p );
		
		regions.push_back(reg);

		regNode = regNode->NextSibling("region");
	}

	regCount = regions.size();

	// чтение параметров о ГРАНИЧНЫХ УСЛОВИЯХ
	node0 = task->FirstChild("boundaries");
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	while (bNode != NULL)
	{
		int edgeType;
		bNode->ToElement()->Attribute("edgeType", &edgeType);
		
		
		HeatBoundary * b = HeatBoundary::create(bNode);
		
		
		
		if (strcmp(str, "BOUND_WALL") == 0) 
		{
			b->parCount = 0;
			b->par = NULL;
			//b->type = Boundary::BOUND_WALL;
		} else
		if (strcmp(str, "BOUND_OUTLET") == 0) 
		{
			b->parCount = 0;
			b->par = NULL;
			//b->type = Boundary::BOUND_OUTLET;
		} else
		if (strcmp(str, "BOUND_INLET") == 0) 
		{
			b->parCount = 4;
			b->par = new double[4];
			//b->type = Boundary::BOUND_INLET;

			node1 = bNode->FirstChild("parameters");
			node1->FirstChild("T")->ToElement()->Attribute("value", &b->par[0]);
			node1->FirstChild("P")->ToElement()->Attribute("value", &b->par[1]);
			node1->FirstChild("Vx")->ToElement()->Attribute("value", &b->par[2]);
			node1->FirstChild("Vy")->ToElement()->Attribute("value", &b->par[3]);
		} else {
			log("ERROR: unsupported boundary condition type '%s'", str);
			EXIT(1);
		}

		boundaries.push_back(b);
		
		
		bNode = bNode->NextSibling("boundCond");
	}
	bCount = boundaries.size();

	node0 = task->FirstChild("mesh");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");
	grid.initFromFiles((char*)fName);


	T		= new double[grid.cCount];
	T_old	= new double[grid.cCount];
	T_int	= new double[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		Region & reg = getRegion(i);
		T[i] = reg.par.T;
	}

	memcpy(T_old, T, grid.cCount*sizeof(double));

	calcTimeStep();
	save(0);
}


void FVM_Heat::calcTimeStep()
{
	double tau = 1.0e+20;
	TAU = _min_(TAU, tau);
	printf("\n\nTime step TAU = %e.\n\n", TAU);
}



void FVM_Heat::run() 
{


	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX) 
	{
		t += TAU; 
		step++;
		memset(T_int, 0, nc*sizeof(double));
		//calcGrad();
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l		= grid.edges[iEdge].l;
			Param pL, pR;
			reconstruct(iEdge, pL, pR);
			double __GAM = 1.4; // TODO: сделать правильное вычисление показателя адиабаты
			calcFlux(fr, fu, fv, fe, pL, pR, n, __GAM);
			
			T_int[c1] += fr*l;
			if (c2 > -1) 
			{
				T_int[c2] -= fr*l;
			}

		}
		memcpy(T, T_old, nc*sizeof(double));
		for (int iCell = 0; iCell < nc; iCell++)
		{
			register double cfl = TAU/grid.cells[iCell].S;
			T[iCell] -= cfl*T_int[iCell];
		}
		memcpy(T_old, T, nc*sizeof(double));
		
		
		if (step % FILE_SAVE_STEP == 0)
		{
			save(step);
		}
		if (step % PRINT_STEP == 0)
		{
			log("step: %d\t\ttime step: %.16f\n", step, t);
		}
	}

}

void FVM_Heat::save(int step)
{
	char fName[50];
	//sprintf(fName, "res_%010d.dat", step);
	//FILE * fp = fopen(fName, "w");
	////fprintf(fp, "# x y ro p u v e M\n");
	//for (int i = 0; i < grid.cCount; i++)
	//{
	//	Param p;
	//	convertConsToPar(i, p);
	//	fprintf(fp, "%25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f %25.16f\n",	grid.cells[i].c.x,
	//																						grid.cells[i].c.y,
	//																						p.r,								//	плотность
	//																						p.p,								//	давление
	//																						p.u,								//	скорость
	//																						p.v,								//
	//																						p.e,								//	внутренняя энергия
	//																						sqrt(_sqr_(p.u)+_sqr_(p.v))/p.cz);	//	число Маха
	//}
	//fclose(fp);
	//// выводим градиенты примитивных переменных
	//printf("File '%s' saved...\n", fName);
	//sprintf(fName, "grad_%010d.dat", step);
	//fp = fopen(fName, "w");
	////fprintf(fp, "# x y ro p u v e M\n");
	//for (int i = 0; i < grid.cCount; i++)
	//{
	//	fprintf(fp, "%25.16f\t%25.16f\t\t%25.16f\t%25.16f\t\t%25.16f\t%25.16f\t\t%25.16f\t%25.16f\n",	primGrad[ID_R][i].x,
	//																									primGrad[ID_R][i].y,
	//																									primGrad[ID_P][i].x,
	//																									primGrad[ID_P][i].y,
	//																									primGrad[ID_U][i].x,
	//																									primGrad[ID_U][i].y,
	//																									primGrad[ID_V][i].x,
	//																									primGrad[ID_V][i].y	);
	//}
	//fclose(fp);
	//printf("File '%s' saved...\n", fName);

	sprintf(fName, "res_%010d.vtk", step);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "GASDIN data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x,  grid.nodes[i].y, 0.0);
		if (i+1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4*grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Density float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.r);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Pressure float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.p);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", p.T);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "SCALARS MachNumber float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f ", (p.u*p.u+p.v*p.v)/p.cz);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}

	fprintf(fp, "VECTORS Velosity float\n");
	for (int i = 0; i < grid.cCount; i++)
	{
		Param p;
		convertConsToPar(i, p);
		fprintf(fp, "%25.16f %25.16f %25.16f ", p.u, p.v, 0.0);
		if (i+1 % 8 == 0 || i+1 == grid.cCount) fprintf(fp, "\n");
	}


	fclose(fp);
	printf("File '%s' saved...\n", fName);


}

void FVM_Heat::calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM)
{
	{	// GODUNOV FLUX
		double RI, EI, PI, UI, VI, WI, UN, UT;
		//double unl = pL.u*n.x+pL.v*n.y;
		//double unr = pR.u*n.x+pR.v*n.y;
		//rim(RI, EI, PI, UN, UI, VI,  pL.r, pL.p, unl, pL.u, pL.v,  pR.r, pR.p, unr, pR.u, pR.v, n, GAM);
			
		double unl = pL.u*n.x+pL.v*n.y;
		double unr = pR.u*n.x+pR.v*n.y;
		double utl = pL.u*n.y-pL.v*n.x;
		double utr = pR.u*n.y-pR.v*n.x;
		rim_orig(RI, EI, PI, UN, UT, WI,  pL.r, pL.p, unl, utl, 0,  pR.r, pR.p, unr, utr, 0, GAM);
		
		UI = UN*n.x+UT*n.y;
		VI = UN*n.y-UT*n.x;

		fr = RI*UN;
		fu = fr*UI+PI*n.x;
		fv = fr*VI+PI*n.y;
		fe = (RI*(EI+0.5*(UI*UI+VI*VI))+PI)*UN;
	}
	//{	// LAX-FRIEDRIX FLUX
	//	double unl = pL.u*n.x+pL.v*n.y;
	//	double unr = pR.u*n.x+pR.v*n.y;
	//	double rol, rul, rvl, rel,  ror, rur, rvr, rer;
	//	double alpha = _max_(fabs(unl)+sqrt(GAM*pL.p/pL.r), fabs(unr)+sqrt(GAM*pR.p/pR.r));
	//	pL.getToCons(rol, rul, rvl, rel);
	//	pR.getToCons(ror, rur, rvr, rer);
	//	double frl = rol*unl;
	//	double frr = ror*unr;
	//	fr = 0.5*(frr+frl								- alpha*(ror-rol));
	//	fu = 0.5*(frr*pR.u+frl*pL.u + (pR.p+pL.p)*n.x	- alpha*(rur-rul));
	//	fv = 0.5*(frr*pR.v+frl*pL.v + (pR.p+pL.p)*n.y	- alpha*(rvr-rvl));
	//	fe = 0.5*((rer+pR.p)*unr + (rel+pL.p)*unl		- alpha*(rer-rel));
	//}
}


void FVM_Heat::reconstruct(int iEdge, Param& pL, Param& pR)
{
	if (grid.edges[iEdge].type == Edge::TYPE_INNER) 
	{
		int c1	= grid.edges[iEdge].c1;
		int c2	= grid.edges[iEdge].c2;
		convertConsToPar(c1, pL);
		convertConsToPar(c2, pR);
	} else {
		int c1	= grid.edges[iEdge].c1;
		convertConsToPar(c1, pL);
		boundaryCond(iEdge, pL, pR);
	}
}


void FVM_Heat::boundaryCond(int iEdge, Param& pL, Param& pR)
{
	int iBound = -1;
	for (int i = 0; i < bCount; i++)
	{
		if (grid.edges[iEdge].type == boundaries[i]->edgeType)
		{
			iBound = i;
			break;
		}
	}
	if (iBound < 0)
	{
		log("ERROR (boundary condition): unknown edge type of edge %d...\n", iEdge);
		EXIT(1);
	}

	boundaries[iBound]->run(pL, pR);

	//Boundary& b = boundaries[iBound];
	//int c1	= grid.edges[iEdge].c1;
	//Material& m = getMaterial(c1);
	//switch (b.type)
	//{
	//case Boundary::BOUND_INLET:
	//	pR.T  = b.par[0];		//!< температура
	//	pR.p  = b.par[1];		//!< давление
	//	pR.u  = b.par[2];		//!< первая компонента вектора скорости
	//	pR.v  = b.par[3];		//!< вторая компонента вектора скорости
	//	
	//	m.URS(pR, 2);
	//	m.URS(pR, 1);
	//	break;
	//
	//case Boundary::BOUND_OUTLET:
	//	pR = pL;
	//	break;
	//
	//case Boundary::BOUND_WALL:
	//	pR = pL;
	//	double Un = pL.u*grid.edges[iEdge].n.x+pL.v*grid.edges[iEdge].n.y;
	//	Vector V;
	//	V.x = grid.edges[iEdge].n.x*Un*2.0; 
	//	V.y = grid.edges[iEdge].n.y*Un*2.0;
	//	pR.u = pL.u-V.x;
	//	pR.v = pL.v-V.y;
	//	break;
	//}
}


void FVM_Heat::done()
{
	delete[] T;
	delete[] T_old;
	delete[] T_int;
}




Region & FVM_Heat::getRegionByCellType(int type)
{
	for (int i = 0; i < regCount; i++)
	{
		if (regions[i].cellType == type) return regions[i];
	}
	log("ERROR: unknown cell type %d...\n", type);
	EXIT(1);
}


Region   &	FVM_Heat::getRegion	(int iCell)
{
	return getRegionByCellType( grid.cells[iCell].type );
}

Material &	FVM_Heat::getMaterial	(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FVM_Heat::convertParToCons(int iCell, Param & par)
{
	T[iCell] = par.T;
}

void FVM_Heat::convertConsToPar(int iCell, Param & par)
{
	par.T = T[iCell];
	//par.u = ru[iCell]/ro[iCell];
	//par.v = rv[iCell]/ro[iCell];
	//par.E = re[iCell]/ro[iCell];
	//par.e = par.E-0.5*(par.u*par.u+par.v*par.v);
	//Material& mat = getMaterial(iCell);
	//mat.URS(par, 0);
}



