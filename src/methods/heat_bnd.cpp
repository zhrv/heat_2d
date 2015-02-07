#include "heat_bnd.h"

const char* HeatBoundary::TYPE_CONST	= "BOUND_CONST";
const char* HeatBoundary::TYPE_FLOW		= "BOUND_FLOW";
const char* HeatBoundary::TYPE_NOFLOW	= "BOUND_NOFLOW";


HeatBoundary* HeatBoundary::create(TiXmlNode* bNode)
{
	HeatBoundary * b = 0;
	TiXmlNode* node = 0;
	TiXmlNode* node1 = 0;

	const char * type = bNode->FirstChild("type")->ToElement()->GetText();
	
	if (strcmp(type, HeatBoundary::TYPE_CONST) == 0) {
		b = new HeatBndConst();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 1;
		b->par = new double[1];
		node = bNode->FirstChild("parameters");
		node1 = node->FirstChild("T");
		if (!node1) throw Exception("Parameter 'T' isn't specified  for BOUND_FLOW.", Exception::TYPE_BOUND_NOPAR); 
		node1->ToElement()->Attribute("value", &b->par[0]);

	}

	if (strcmp(type, HeatBoundary::TYPE_NOFLOW) == 0) {
		b = new HeatBndNoFlow();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 0;
		b->par = NULL;
	}

	if (strcmp(type, HeatBoundary::TYPE_FLOW) == 0) {
		b = new HeatBndFlow();
		bNode->ToElement()->Attribute("edgeType", &b->edgeType);
		b->parCount = 2;
		b->par = new double[2];
		node = bNode->FirstChild("parameters");
		node1 = node->FirstChild("factor");
		if (!node1) throw Exception("Parameter 'factor' isn't specified  for BOUND_FLOW.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[0]);   // factor*(dT/dX)=flowValue
		
		node1 = node->FirstChild("flowValue");
		if (!node1) throw Exception("Parameter 'flowValue' isn't specified  for BOUND_FLOW.", Exception::TYPE_BOUND_NOPAR);
		node1->ToElement()->Attribute("value", &b->par[1]);

	}

	if (!b) {
		throw Exception("Unknown bountary type '%s' specified.", Exception::TYPE_BOUND_UNKNOWN);
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

