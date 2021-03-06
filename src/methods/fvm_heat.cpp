#include "fvm_heat.h"
#include "tinyxml.h"
#include <string>
#include "global.h"
#include "MeshReader.h"

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
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_SAVE_STEP);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &PRINT_STEP);

	// ������ ���������� � ����������
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
		node1->FirstChild("Kx")->ToElement()->Attribute("value", &mat.Kx);
		node1->FirstChild("Ky")->ToElement()->Attribute("value", &mat.Ky);
		materials.push_back(mat);

		matNode = matNode->NextSibling("material");
	}
	matCount = materials.size();

	// ������ ���������� � ��������
	node0 = task->FirstChild("regions");
	TiXmlNode* regNode = node0->FirstChild("region");
	while (regNode != NULL)
	{
		Region reg;
		regNode->ToElement()->Attribute("id", &reg.id);
		regNode->FirstChild("material")->ToElement()->Attribute("id", &reg.matId);
		regNode->FirstChild("cell")->ToElement()->Attribute("type", &reg.cellType);
		
		reg.name = regNode->FirstChild("name")->ToElement()->GetText();

		node1 = regNode->FirstChild("parameters");
		node1->FirstChild("T")->ToElement()->Attribute("value", &reg.par.T);
		node1->FirstChild("P")->ToElement()->Attribute("value", &reg.par.p);

		node1->FirstChild("Vx")->ToElement()->Attribute("value", &reg.par.u);
		node1->FirstChild("Vy")->ToElement()->Attribute("value", &reg.par.v);

		regions.push_back(reg);

		regNode = regNode->NextSibling("region");
	}

	regCount = regions.size();

	/* ������ ���������� �� */
	node0 = task->FirstChild("boundaries");
	TiXmlNode* bNode = node0->FirstChild("boundCond");
	while (bNode != NULL)
	{
		int edgeType;
		bNode->ToElement()->Attribute("edgeType", &edgeType);
		
		HeatBoundary * b;

		try {
			b = HeatBoundary::create(bNode);
		}
		catch (Exception e) {
			log("ERROR: %s\n", e.getMessage());
			exit(e.getType());
		}
		
		boundaries.push_back(b);
		
		bNode = bNode->NextSibling("boundCond");
	}
	
	bCount = boundaries.size();

	/* ������ ������ �����. */
	node0 = task->FirstChild("mesh");
	const char* fName = node0->FirstChild("name")->ToElement()->Attribute("value");
	const char* tName = node0->FirstChild("filesType")->ToElement()->Attribute("value");
	MeshReader* mr = MeshReader::create(MeshReader::getType((char*)tName), (char*)fName);
	mr->read(&grid);

	/* ����������� �� ��� ������ ������. */
	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		Edge & e = grid.edges[iEdge];
		if (e.type == Edge::TYPE_INNER) {
			e.bnd = NULL;
			continue;
		}
		if (e.type == Edge::TYPE_NAMED) {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (strcmp(e.typeName, boundaries[i]->name) == 0)
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

			e.bnd = boundaries[iBound];
		}
		else {
			int iBound = -1;
			for (int i = 0; i < bCount; i++)
			{
				if (e.type == boundaries[i]->edgeType)
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

			e.bnd = boundaries[iBound];
		}
	}

	K		= new Vector[grid.cCount];
	T		= new double[grid.cCount];
	T_old	= new double[grid.cCount];
	T_int	= new double[grid.cCount];
	gradT	= new Vector[grid.cCount];
	V		= new Vector[grid.cCount];

	for (int i = 0; i < grid.cCount; i++)
	{
		if (grid.cells[i].type == -1) {
			Region & reg = getRegion(grid.cells[i].typeName);
			T[i] = reg.par.T;
			V[i].x = reg.par.u;
			V[i].y = reg.par.v;
			K[i].x = materials[reg.matId].Kx;
			K[i].y = materials[reg.matId].Ky;
		}
		else {
			Region & reg = getRegion(i);
			T[i] = reg.par.T;
			V[i].x = reg.par.u;
			V[i].y = reg.par.v;
			K[i].x = materials[reg.matId].Kx;
			K[i].y = materials[reg.matId].Ky;
		}
		
	}

	memcpy(T_old, T, grid.cCount*sizeof(double));

	calcTimeStep();
	save(0);
}


void FVM_Heat::calcTimeStep()
{
	double tau = 1.0e+20;
	for (int i = 0; i < grid.cCount; i++) {
		double t = CFL*grid.cells[i].S;
		if (t < tau) tau = t;
	}
	TAU = _min_(TAU, tau);
	log("\n\nTime step TAU = %e.\n\n", TAU);
}

void FVM_Heat::done()
{
	delete[] K;
	delete[] T;
	delete[] T_old;
	delete[] T_int;
	delete[] gradT;
	delete[] V;
}




void FVM_Heat::run() 
{


	int nc = grid.cCount;
	int ne = grid.eCount;

	double			t		= 0.0;
	unsigned int	step	= 0;
	while (t < TMAX && step < STEP_MAX) 
	{
		t += TAU; 
		step++;
		calcGrad();
		memset(T_int, 0, nc*sizeof(double));
		for (int iEdge = 0; iEdge < ne; iEdge++)
		{
			Edge &edge = grid.edges[iEdge];
			double fr, fu, fv, fe;
			int c1	= grid.edges[iEdge].c1;
			int c2	= grid.edges[iEdge].c2;
			Vector n	= grid.edges[iEdge].n;
			double l		= grid.edges[iEdge].l;
			Param pL, pR;
			double q, h;
			if (edge.type == Edge::TYPE_INNER)
			{
				convertConsToPar(c1, pL);
				convertConsToPar(c2, pR);
				q = scalar_prod(gradT[c1], n) + scalar_prod(gradT[c2], n);
				q *= 0.5;
				
				Vector Vi;
				Vi.x = 0.5*(pL.u*pL.T + pR.u*pR.T);
				Vi.y = 0.5*(pL.v*pL.T + pR.v*pR.T);
				q -= scalar_prod(Vi, n);
			}
			else {
				int c1 = grid.edges[iEdge].c1;
				convertConsToPar(c1, pL);
				boundaryCond(iEdge, pL, pR);
				//h = edge.cnl1 + edge.cnl1;
				//q = (pR.T - pL.T) / h;
				q = scalar_prod(gradT[c1], n);
				
				Vector Vi;
				Vi.x = 0.5*(pL.u*pL.T + pR.u*pR.T);
				Vi.y = 0.5*(pL.v*pL.T + pR.v*pR.T);
				q -= scalar_prod(Vi, n);
			}
			T_int[c1] += q*l;
			if (c2 > -1) 
			{
				T_int[c2] -= q*l;
			}

		}
		//memcpy(T, T_old, nc*sizeof(double));
		for (int iCell = 0; iCell < nc; iCell++)
		{
			register double cfl = TAU/grid.cells[iCell].S;
			T[iCell] += cfl*T_int[iCell];
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


void FVM_Heat::calcGrad()
{
	for (int i = 0; i < grid.cCount; i++) {
		gradT[i].x = 0.0;
		gradT[i].y = 0.0;
	}

	for (int iEdge = 0; iEdge < grid.eCount; iEdge++) {
		Edge &edge = grid.edges[iEdge];
		double fr, fu, fv, fe;
		int c1 = grid.edges[iEdge].c1;
		int c2 = grid.edges[iEdge].c2;
		Vector n = grid.edges[iEdge].n;
		double l = grid.edges[iEdge].l;
		Param pL, pR;
		double q, h;
		if (edge.type == Edge::TYPE_INNER)
		{
			convertConsToPar(c1, pL);
			convertConsToPar(c2, pR);
			q = 0.5*(pL.T+pR.T);
		}
		else {
			convertConsToPar(c1, pL);
			boundaryCond(iEdge, pL, pR);
			q = 0.5*(pL.T + pR.T);
		}
		gradT[c1].x += q*l*n.x*K[c1].x;
		gradT[c1].y += q*l*n.y*K[c1].y;
		if (c2 > -1)
		{
			gradT[c2].x -= q*l*n.x*K[c2].x;
			gradT[c2].y -= q*l*n.y*K[c2].y;
		}
	}

	for (int i = 0; i < grid.cCount; i++) {
		gradT[i].x /= grid.cells[i].S;
		gradT[i].y /= grid.cells[i].S;
	}
}

void FVM_Heat::boundaryCond(int iEdge, Param& pL, Param& pR)
{
	if (grid.edges[iEdge].bnd) {
		grid.edges[iEdge].bnd->run(pL, pR);
	}
}


void FVM_Heat::save(int step)
{
	char fName[50];

	sprintf(fName, "res_%010d.vtk", step);
	FILE * fp = fopen(fName, "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n");
	fprintf(fp, "HEAT_2D data file\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(fp, "POINTS %d float\n", grid.nCount);
	for (int i = 0; i < grid.nCount; i++)
	{
		fprintf(fp, "%f %f %f  ", grid.nodes[i].x, grid.nodes[i].y, 0.0);
		if (i + 1 % 8 == 0) fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
	fprintf(fp, "CELLS %d %d\n", grid.cCount, 4 * grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "3 %d %d %d\n", grid.cells[i].nodesInd[0], grid.cells[i].nodesInd[1], grid.cells[i].nodesInd[2]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "CELL_TYPES %d\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++) fprintf(fp, "5\n");
	fprintf(fp, "\n");

	fprintf(fp, "CELL_DATA %d\nSCALARS Temperature float 1\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f ", T[i]);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}


	fprintf(fp, "SCALARS Heat_K float 2\nLOOKUP_TABLE default\n", grid.cCount);
	for (int i = 0; i < grid.cCount; i++)
	{
		fprintf(fp, "%25.16f %25.16f ", K[i].x, K[i].y);
		if (i + 1 % 8 == 0 || i + 1 == grid.cCount) fprintf(fp, "\n");
	}



	fclose(fp);
	printf("File '%s' saved...\n", fName);


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

Region & FVM_Heat::getRegionByName(char* name)
{
	for (int i = 0; i < regCount; i++)
	{
		if (strcmp(regions[i].name, name) == 0) return regions[i];
	}
	log("ERROR: unknown cell name '%d'...\n", name);
	EXIT(1);
}


Region   &	FVM_Heat::getRegion(int iCell)
{
	return getRegionByCellType(grid.cells[iCell].type);
}

Region   &	FVM_Heat::getRegion(char * name)
{
	return getRegionByName(name);
}

Material &	FVM_Heat::getMaterial(int iCell)
{
	Region & reg = getRegion(iCell);
	return materials[reg.matId];
}


void FVM_Heat::convertParToCons(int iCell, Param & par)
{
	T[iCell] = par.T;
	V[iCell].x = par.u;
	V[iCell].y = par.v;
}

void FVM_Heat::convertConsToPar(int iCell, Param & par)
{
	par.T = T[iCell];
	par.u = V[iCell].x;
	par.v = V[iCell].y;
}



