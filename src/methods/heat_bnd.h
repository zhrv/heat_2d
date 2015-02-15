#ifndef _HEAT_BND_H_
#define _HEAT_BND_H_

#include <cstdlib>
#include "global.h"
#include "tinyxml.h"

struct HeatBoundary
{
	static HeatBoundary* create(TiXmlNode* bNode);
	virtual void run(Param& pL, Param& pR) = 0;

	HeatBoundary() : par(NULL) {}
	~HeatBoundary() { if (par) delete[] par; par = NULL; }
	
	
	int			edgeType;
	char		name[64];
	double *	par;
	int			parCount;

	static const char* TYPE_CONST;
	static const char* TYPE_FLOW;
	static const char* TYPE_NOFLOW;

};


struct HeatBndConst : public HeatBoundary
{
	virtual void run(Param& pL, Param& pR);
};

struct HeatBndFlow : public HeatBoundary
{
	virtual void run(Param& pL, Param& pR);
};

struct HeatBndNoFlow : public HeatBoundary
{
	virtual void run(Param& pL, Param& pR);
};




#endif