#ifndef _FVM_TVD_H_
#define _FVM_TVD_H_

#include <method.h>
#include <vector>
#include "heat_bnd.h"

typedef std::vector<Material> Materials;
typedef std::vector<Region> Regions;
typedef std::vector<HeatBoundary*> HeatBoundaries;


class FVM_Heat: public Method
{
public:
	virtual void init(char * xmlFileName);
	virtual void run();
	virtual void done();
protected:
	Region & getRegionByCellType(int type);
	Region & getRegionByName(char* name);

	Region   &	getRegion(int iCell);
	Region   &	getRegion(char* name);
	Material &	getMaterial(int iCell);
	
	/**
	 *	ѕреобразование примитивных переменных в консервативные
	 */
	void convertParToCons(int iCell, Param & par);
	
	/**
	 *	ѕреобразование консервативных переменных в примитивные
	 */
	void convertConsToPar(int iCell, Param & par);
	
	/**
	 *	¬ычисление параметров справа и слева от границы €чейки
	 */
	void reconstruct(int iFace, Param& pL, Param& pR);

	/**
	 *	¬ычисление параметров с внешней стороны от границы €чейки согласно граничным услови€м
	 */
	void boundaryCond(int iFace, Param& pL, Param& pR);

	/**
	 *	¬ычисление шага по времени по значению CFL, если значение TAU из XML 
	 *	меньше вычисленного, то используетс€ значение, заданное в XML
	 */
	void calcTimeStep();

	/**
	 *	«апись значений газодинамических параметров в файл
	 */
	void save(int);

	/**
	*	¬ычисление численного потока
	*/
	void calcFlux(double& fr, double& fu, double& fv, double& fe, Param pL, Param pR, Vector n, double GAM);

	/**
	*	¬ычисление градиентов
	*/
	void calcGrad();

private:
	double TMAX;
	double TAU;
	double CFL;
	int FILE_SAVE_STEP;
	int PRINT_STEP;
	int STEP_MAX;

	int				matCount;
	int				regCount;
	int				bCount;
	Materials		materials;
	Regions			regions;
	HeatBoundaries		boundaries;

	Vector * K;
	double * T;
	double * T_old;
	double * T_int;
	Vector * gradT;
	Vector * V;

};

#endif