#include <stdlib.h>
#include <stdio.h>
#include <solver.h>
#include "global.h"
#include <cfloat>

int main(int argc, char** argv)
{
#ifdef _DEBUG
	_controlfp(~(_MCW_EM & (~_EM_INEXACT) & (~_EM_UNDERFLOW)), _MCW_EM);
#endif

	hLog = fopen("task.log", "w"); // ��������� ���� ��� ������ ����; ������ printf(...) ���������� ������������ log(...)
	
	Method * method = Solver::initMethod( "task.xml" ); 
	
	Solver :: runMethod		( method );
	Solver :: destroyMethod	( method );

	fclose(hLog);

	return 0;
}