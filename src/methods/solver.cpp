#include <solver.h>
#include <method.h>
#include <fvm_heat.h>
#include <decomp.h>
#include "tinyxml.h"
#include <string.h>



Method* Solver::initMethod(char* fileName)
{
	Method * m;

	TiXmlDocument doc( fileName );
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
	//char methodName[50];
	const char * methodName = task->ToElement()->Attribute("method");

	if (strcmp("Decomp", methodName) == 0)
	{
		m = new Decomp();
	}
	else if (strcmp("FVM_HEAT", methodName) == 0) 
	{
		m = new FVM_Heat();
	}
	else
	{
		log("ERROR: unknown method '%s'.\n", methodName);
		exit(-1);
	}

	m->init(fileName);
	return m;
}

void Solver::runMethod(Method* m)
{
	m->run();
}

void Solver::destroyMethod(Method* m)
{
	m->done();
	delete m;
}