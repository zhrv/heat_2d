#pragma once
#include "MeshReader.h"
class MeshReaderSalome :
	public MeshReader
{
private:
	char* fileName;
public:
	MeshReaderSalome(char* fName) : fileName(fName) {}
	~MeshReaderSalome() { delete[] fileName; }

	virtual void read(Grid*);
};

