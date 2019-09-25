#ifndef _MYWORLD_
#define _MYWORLD_

#include <vector>
//#include <Eigen\IterativeLinearSolvers>

#define IX(i, j) ((i)+(getNumCells()+2)*(j))
#define SWAPPING_1D(x0,x) {double *tmp=x0;x0=x;x=tmp;}
//#define SWAPPING_2D(x0,x) {double **tmp=x0;x0=x;x=tmp;}
#define SWAPPING_2D(x0,x) {double tmp[3]={1,2,3};tmp=x0;x0=x;x=tmp;}
#define SIZE 4
#define SIZE_I 4
#define SIZE_J 4


class MyWorld {
 public:
    MyWorld();

    virtual ~MyWorld();

    void initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef);

    double getTimeStep();
    int getNumCells() { return mNumCells;}
	int getNumCells_I() { return SIZE_I; }
	int getNumCells_J() { return SIZE_J; }

    //double getDensity(int _index) { return mDensity[_index]; }
	/*double getDensity_R(int _index) { return mDensity_R[_index]; }
	double getDensity_G(int _index) { return mDensity_G[_index]; }
	double getDensity_B(int _index) { return mDensity_B[_index]; }*/
	double getDensity_R(int _indexI, int _indexJ) { return mDensity_R[_indexI][_indexJ]; }
	double getDensity_G(int _indexI, int _indexJ) { return mDensity_G[_indexI][_indexJ]; }
	double getDensity_B(int _indexI, int _indexJ) { return mDensity_B[_indexI][_indexJ]; }

    /*double getVelocityU(int _index) { return mU[_index]; }
    double getVelocityV(int _index) { return mV[_index]; }*/
	double getVelocityU(int _indexI, int _indexJ) { return mU[_indexI][_indexJ]; }
	double getVelocityV(int _indexI, int _indexJ) { return mV[_indexI][_indexJ]; }

	double getDivergence(int _indexI, int _indexJ) { return mV[_indexI][_indexJ]; }
	double getAverageDivergence() { return averageDivergence; }
	double getAverageAbsDivergence() { return averageAbsDivergence; }

    //void setDensity(int _i, int _j, double _source) { mDensity[IX(_i, _j)] += mTimeStep * _source; }
	//void setDensity_R(int _i, int _j, double _source) { mDensity_R[IX(_i, _j)] += mTimeStep * _source; }
	//void setDensity_G(int _i, int _j, double _source) { mDensity_G[IX(_i, _j)] += mTimeStep * _source; }
	//void setDensity_B(int _i, int _j, double _source) { mDensity_B[IX(_i, _j)] += mTimeStep * _source; }
	void setDensity_R(int _i, int _j, double _source) { mDensity_R[_i][_j] += mTimeStep * _source; }
	void setDensity_G(int _i, int _j, double _source) { mDensity_G[_i][_j] += mTimeStep * _source; }
	void setDensity_B(int _i, int _j, double _source) { mDensity_B[_i][_j] += mTimeStep * _source; }

	/*void setU(int _i, int _j, double _force) { mU[IX(_i, _j)] += mTimeStep * _force; }
    void setV(int _i, int _j, double _force) { mV[IX(_i, _j)] += mTimeStep * _force; }*/
	void setU(int _i, int _j, double _force) { mU[_i][_j] += mTimeStep * _force; }
	void setV(int _i, int _j, double _force) { mV[_i][_j] += mTimeStep * _force; }
    
    void simulate();
    
 protected:
	void densityStep_1D(double *_x, double *_x0);
	/*void densityStep_2D(double **_x, double **_x0);*/
	void densityStep_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2]);

    void velocityStep_1D(double *_u, double *_v, double *_u0, double *_v0);
	void velocityStep_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]);

    void diffuseDensity_1D(double *_x, double *_x0);
	//void diffuseDensity_2D(double **_x, double **_x0);
	void diffuseDensity_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2]);

	void diffuseVelocity_1D(double *_u, double *_v, double *_u0, double *_v0);
	void diffuseVelocity_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]);

    void advectDensity_1D(double *_d, double *_d0, double *_u, double *_v);
	//void advectDensity_2D(double **_d, double **_d0, double *_u, double *_v);
	//void advectDensity_2D(double _d[][SIZE + 2], double _d0[][SIZE + 2], double *_u, double *_v);
	void advectDensity_2D(double _d[][SIZE_J + 2], double _d0[][SIZE_J + 2], double _u[][SIZE_J + 2], double _v[][SIZE_J + 2]);

    void advectVelocity_1D(double *_u, double *_v, double *_u0, double *_v0);
	void advectVelocity_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]);

    void project_1D(double *_u, double *_v, double *_u0, double *_v0);
	void project_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]);

    /*void externalForces_1D();*/
	void externalForces_2D();

	void calcDivergence(double _d[][SIZE_J + 2]);

    void linearSolve_1D(double *_x, double *_x0, double _a, double _c);
	//void linearSolve_2D(double **_x, double **_x0, double _a, double _c);
	void linearSolve_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2], double _a, double _c);

    void setBoundary_1D(double *_x);
	/*void setBoundary_2D(double **_x);*/
	void setBoundary_2D(double _x[][SIZE_J + 2]);

    void setVelocityBoundary_1D(double *_u, double *_v);
	void setVelocityBoundary_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2]);

	void initializeCG();
	void project_CG(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _d[][SIZE_J + 2]);

	//void swap_2D(double _x[][SIZE + 2], double _x0[][SIZE + 2]);

    int mNumCells;
	int sizeI;
	int sizeJ;
    double mTimeStep;
    double mDiffusionCoef;
    double mViscosityCoef;

    /*double *mU;
    double *mV;
    double *mPreU;
    double *mPreV;*/

	double mU[SIZE_I + 2][SIZE_J + 2];
	double mV[SIZE_I + 2][SIZE_J + 2];
	double mPreU[SIZE_I + 2][SIZE_J + 2];
	double mPreV[SIZE_I + 2][SIZE_J + 2];

    //double *mDensity;
    //double *mPreDensity;

	/*double *mDensity_R;
	double *mPreDensity_R;
	double *mDensity_G;
	double *mPreDensity_G;
	double *mDensity_B;
	double *mPreDensity_B;*/

	//double myArray [2][3];

	double mDensity_R[SIZE_I + 2][SIZE_J + 2];
	double mPreDensity_R[SIZE_I + 2][SIZE_J + 2];
	double mDensity_G[SIZE_I + 2][SIZE_J + 2];
	double mPreDensity_G[SIZE_I + 2][SIZE_J + 2];
	double mDensity_B[SIZE_I + 2][SIZE_J + 2];
	double mPreDensity_B[SIZE_I + 2][SIZE_J + 2];

	double mDivergenceArray[SIZE_I + 2][SIZE_J + 2];
	double averageDivergence;
	double averageAbsDivergence;

	double temp2DArray[SIZE + 2][SIZE + 2];

	//EIGEN_CONJUGATE_GRADIENT_H<

	/*double **mDensity_R;
	double **mPreDensity_R;
	double **mDensity_G;
	double **mPreDensity_G;
	double **mDensity_B;
	double **mPreDensity_B;*/
};

#endif
