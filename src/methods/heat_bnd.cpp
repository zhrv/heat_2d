#include "heat_bnd.h"

const char* HeatBoundary::TYPE_CONST	= "BOUND_CONST";
const char* HeatBoundary::TYPE_FLOW		= "BOUND_FLOW";
const char* HeatBoundary::TYPE_NOFLOW	= "BOUND_NOFLOW";


HeatBoundary* HeatBoundary::create(TiXmlNode* bNode)
{
	HeatBoundary * b;
	TiXmlNode* node = 0;

	const char * type = bNode->FirstChild("type")->ToElement()->GetText();
	
	if (strcmp(type, HeatBoundary::TYPE_CONST) == 0) {
		b = new HeatBndConst();
		b->parCount = 1;
		b->par = new double[1];
		node = bNode->FirstChild("parameters");
		node->FirstChild("T")->ToElement()->Attribute("value", &b->par[0]);

	}

	if (strcmp(type, HeatBoundary::TYPE_NOFLOW) == 0) {
		b = new HeatBndNoFlow();
		b->parCount = 0;
		b->par = NULL;
	}

	if (strcmp(type, HeatBoundary::TYPE_FLOW) == 0) {
		b = new HeatBndConst();
		b->parCount = 2;
		b->par = new double[2];
		node = bNode->FirstChild("parameters");
		node->FirstChild("factor")->ToElement()->Attribute("value", &b->par[0]);   // factor*(dT/dX)=flowValue
		node->FirstChild("flowValue")->ToElement()->Attribute("value", &b->par[1]);

	}

	return b;
}


void HeatBndConst::run(Param& pL, Param& pR)
{
	pR.T = par[0];
}

void HeatBndFlow::run(Param& pL, Param& pR) // @todo исправить - необходимо учитывать шаг сетки
{
	pR.T = pL.T+par[1]/par[0];
}

void HeatBndNoFlow::run(Param& pL, Param& pR)
{
	pR.T = pL.T;
}

