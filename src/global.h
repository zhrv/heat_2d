#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
	#include <direct.h>
	#define MK_DIR(name) _mkdir("mesh");
#else
	#include <sys/stat.h>
	#define MK_DIR(name) mkdir(name, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif



//! ������������� ������� ����������
const double gR = 8.314472;	// �2 �� �-2 �-1 ����-1

/**
 *	����� �� ���������
 */
struct Point
{
	double x;
	double y;
};


/**
 *	��������� ������
 */
typedef Point Vector;

/**
 *	������ ����������� n
 */
struct VECTOR
{
	VECTOR(int an = 0): n(an) { if (n) { elem = new double[n]; memset(elem, 0, sizeof(double)*n); } else elem = NULL; }
	VECTOR(const VECTOR& v): n(v.n) { elem = new double[n]; memcpy(elem, v.elem, n*sizeof(double)); }
	~VECTOR() { if (elem) delete[] elem; n = 0; }

	void init(int an) { n = an; if (elem) delete[] elem; elem = new double[n]; }

	double&		operator []		(int i)					{ return elem[i]; }
	
	double&		operator ()		(int i)					{ return elem[i]; }
	
	void		operator =		(const VECTOR& v)		{	
		if (elem) delete[] elem; 
		n = v.n; 
		elem = new double[n]; 
		memcpy(elem, v.elem, n*sizeof(double)); 
	}

	void		operator =		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] = x; }

	void		operator +=		(const VECTOR& v)		{ for (int i = 0; i < n; i++) elem[i] += v.elem[i]; }	

	void		operator +=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] += x;	}	

	void		operator -=		(const VECTOR& v)		{ for (int i = 0; i < n; i++) elem[i] -= v.elem[i];	}		
	
	void		operator -=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] -= x; }

	void		operator *=		(const double& x)		{ for (int i = 0; i < n; i++) elem[i] *= x;	}		

	static double SCALAR_PROD(const VECTOR& v1, const VECTOR& v2) {
		double s = 0;
		for (int i = 0; i < v1.n; i++)
		{
			s += v1.elem[i]*v2.elem[i];
		}
		return s;
	}

	double *	elem;
	int			n;
};



/**
 * ������� ����������� nxn
 */
struct MATRIX
{
	MATRIX(int an): n(an), count(an*an) { elem = new double[n*n]; memset(elem, 0, sizeof(double)*count); }
	MATRIX(const MATRIX& m): n(m.n), count(m.n*m.n) { if (elem) delete[] elem; elem = new double[count]; memcpy(elem, m.elem, count*sizeof(double)); }
	~MATRIX() { if (elem) delete[] elem; n = 0; }

	void operator = (const MATRIX& m) {
		if (elem) delete[] elem; 
		n = m.n;
		elem = new double[count]; 
		memcpy(elem, m.elem, count*sizeof(double));
	}

	MATRIX&		operator  =		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] = x;	}	

	MATRIX&		operator +=		(const MATRIX& v)		{ for (int i = 0; i < count; i++) elem[i] += v.elem[i]; }	

	MATRIX&		operator +=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] += x;	}	

	MATRIX&		operator -=		(const MATRIX& v)		{ for (int i = 0; i < count; i++) elem[i] -= v.elem[i];	}		
	
	MATRIX&		operator -=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] -= x; }

	MATRIX&		operator *=		(const double& x)		{ for (int i = 0; i < count; i++) elem[i] *= x; }		


	double *	elem;
	int			n;
	int			count;
};

/**
 *	���������������� ���������
 */
struct Param
{
	double r;		//!< ���������
	double p;		//!< ��������
	double e;		//!< ���������� �������
	double E;		//!< ������ �������
	double u;		//!< ������ ���������� ������� ��������
	double v;		//!< ������ ���������� ������� ��������
	double cz;		//!< �������� �����
	double T;		//!< �����������
	double ML;		//!< ������������ ��������
	
};

/**
 *	��������� ������� ��������� ������
 */
struct Region
{
	int		id;
	int		matId;
	int		cellType;		//!< ��� ������
	Param	par;			//!< ��������� �������
	const char*	name;		//!< ��� �������
};

struct Material
{
	const char*	name;
	int			id;
	
	double		M;			//!< �������� �����
	double		Cp;			//!< ������������ ��� ���������� ��������
	double		ML;			//!< ������������ ��������
	double		Kx, Ky;		//!< ����������� ����������������

	void URS(Param &par, int flag);
};

struct Boundary
{
	Boundary(): par(NULL) {}
	~Boundary() { if (par) delete[] par; par = NULL; }
	int			edgeType;
	int			type;
	double *	par;
	int			parCount;

	static const int BOUND_INLET	= 1;
	static const int BOUND_OUTLET	= 2;
	static const int BOUND_WALL		= 3;
};


class Exception
{
public:
	static const int TYPE_BOUND_NOPAR = 101;
	static const int TYPE_BOUND_UNKNOWN = 102;

	static const int TYPE_MESH_WRONG_NAME = 201;
	static const int TYPE_MESH_UNV_UNKNOWN_ELEMENT = 221;
	static const int TYPE_MESH_UNV_NOT_DEFINED_BNT_EDGE = 222;

	static const int FILE_OPENING_ERROR = 301;

	Exception(char* msg, int t) : message(msg), type(t) {}
	char* getMessage() { return message; }
	int getType() { return type; }

private:
	char* message;
	int type;
};


extern FILE * hLog;

inline double scalar_prod(Vector a, Vector b) { return a.x*b.x+a.y*b.y; }

extern void log(char * format, ...);
extern void EXIT(int err);

/**
 *	������� ������ � ������� ������������� �������
 *	(�) ��� ��. �.�.������� ���, ������, �������, ���������
 *
 *	c==========================================================
 *	C    Nikichine
 *	C    module for tube.for /Zmitrenko/
 *	c==========================================================
 */
extern void rim(double& RI, double& EI, double& PI, double& UN, double& UI, double& VI,
         double RB, double PB, double UNB, double UB, double VB,
         double RE, double PE, double UNE, double UE, double VE, Vector n, double GAM);


extern void rim_orig(double& RI, double& EI, double& PI, double& UI, double& VI, double& WI,
         double RB, double PB, double UB, double VB, double WB,
         double RE, double PE, double UE, double VE, double WE, double GAM);

#endif