#include "MyWorld.h"
#include <string>
#include <iostream>
#include <stdlib.h>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

ConjugateGradient<SparseMatrix<double>, Lower | Upper> cg;
SparseMatrix<double> matA(SIZE_I * SIZE_J, SIZE_I * SIZE_J);

int frameCount = 0;

MyWorld::MyWorld() {
}

MyWorld::~MyWorld() {
}

void MyWorld::initialize(int _numCells, double _timeStep, double _diffCoef, double _viscCoef) {

	initializeCG();
	
    mNumCells = SIZE;
    mTimeStep = _timeStep;
    mDiffusionCoef = _diffCoef;
    mViscosityCoef = _viscCoef;

    int size = (SIZE + 2) * (SIZE + 2);
	sizeI = SIZE_I + 2;
	sizeJ = SIZE_J + 2;

    // Allocate memory for velocity and density fields
    /*mU = new double[size];
    mV = new double[size];
    mPreU = new double[size];
    mPreV = new double[size];*/
    //mDensity = new double[size];
    //mPreDensity = new double[size];

	/*mDensity_R = new double[size];
	mPreDensity_R = new double[size];
	mDensity_G = new double[size];
	mPreDensity_G = new double[size];
	mDensity_B = new double[size];
	mPreDensity_B = new double[size];*/


	//// Set up density arrays as 2D arrays
	//mDensity_R = new double*[sizeI];
	//mPreDensity_R = new double*[sizeI];
	//mDensity_G = new double*[sizeI];
	//mPreDensity_G = new double*[sizeI];
	//mDensity_B = new double*[sizeI];
	//mPreDensity_B = new double*[sizeI];
	//for (int i = 0; i < sizeI; i++) {
	//	mDensity_R[i] = new double[sizeJ];
	//	mPreDensity_R[i] = new double[sizeJ];
	//	mDensity_G[i] = new double[sizeJ];
	//	mPreDensity_G[i] = new double[sizeJ];
	//	mDensity_B[i] = new double[sizeJ];
	//	mPreDensity_B[i] = new double[sizeJ];
	//}
    
  //  for (int i = 0; i < size; i++) {
  //      mU[i] = mPreU[i] = 0.0;
  //      mV[i] = mPreV[i] = 0.0;
  //      //mDensity[i] = mPreDensity[i] = 0.0;
		/*mDensity_R[i] = mPreDensity_R[i] = 0.0;
		mDensity_G[i] = mPreDensity_G[i] = 0.0;
		mDensity_B[i] = mPreDensity_B[i] = 0.0;*/
  //  }  

	// Zero out arrays (2D arrays)
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			mDensity_R[i][j] = mPreDensity_R[i][j] = 0.0;
			mDensity_G[i][j] = mPreDensity_G[i][j] = 0.0;
			mDensity_B[i][j] = mPreDensity_B[i][j] = 0.0;
			mU[i][j] = mPreU[i][j] = 0.0;
			mV[i][j] = mPreV[i][j] = 0.0;
			mDivergenceArray[i][j] = 0.0;
			averageDivergence = 0.0;
			//temp2DArray[i][j] = 0.0;
		}
	}
	mU[2][2] = 1.0;

}

// Initializes the A matrix for the pressure solve
void MyWorld::initializeCG() {
	std::vector<Triplet<double>> tripletList;
	int numOfEntries = 0;
	numOfEntries += 3 * 4;
	numOfEntries += 4 * (SIZE_J - 2) * 2;
	numOfEntries += 4 * 2 * (SIZE_J - 2);
	numOfEntries += 5 * (SIZE_I - 2) * (SIZE_J - 2);
	tripletList.reserve(numOfEntries);
	for (int i = 0; i < SIZE_I * SIZE_J; i++) {
		if (i == 0) {
			tripletList.push_back(Triplet<double>(i, i, -2));
			tripletList.push_back(Triplet<double>(i, i + 1, 1));
			tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
		}
		else if (i > 0 && i < SIZE_I - 1) {
			tripletList.push_back(Triplet<double>(i, i, -3));
			tripletList.push_back(Triplet<double>(i, i - 1, 1));
			tripletList.push_back(Triplet<double>(i, i + 1, 1));
			tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
		}
		else if (i == SIZE_I - 1) {
			tripletList.push_back(Triplet<double>(i, i, -2));
			tripletList.push_back(Triplet<double>(i, i - 1, 1));
			tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
		}
		else if (i == (SIZE_I - 1) * SIZE_J) {
			tripletList.push_back(Triplet<double>(i, i, -2));
			tripletList.push_back(Triplet<double>(i, i + 1, 1));
			tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
		}
		else if (i > (SIZE_I - 1) * SIZE_J && i < (SIZE_I * SIZE_J) - 1) {
			tripletList.push_back(Triplet<double>(i, i, -3));
			tripletList.push_back(Triplet<double>(i, i - 1, 1));
			tripletList.push_back(Triplet<double>(i, i + 1, 1));
			tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
		}
		else if (i == (SIZE_I * SIZE_J) - 1) {
			tripletList.push_back(Triplet<double>(i, i, -2));
			tripletList.push_back(Triplet<double>(i, i - 1, 1));
			tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
		}
		else {
			if (i % SIZE_J == 0) {
				tripletList.push_back(Triplet<double>(i, i, -3));
				tripletList.push_back(Triplet<double>(i, i + 1, 1));
				tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
				tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
			}
			else if (i % SIZE_J == SIZE_J - 1) {
				tripletList.push_back(Triplet<double>(i, i, -3));
				tripletList.push_back(Triplet<double>(i, i - 1, 1));
				tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
				tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
			}
			else {
				tripletList.push_back(Triplet<double>(i, i, -4));
				tripletList.push_back(Triplet<double>(i, i - 1, 1));
				tripletList.push_back(Triplet<double>(i, i + 1, 1));
				tripletList.push_back(Triplet<double>(i, i - SIZE_J, 1));
				tripletList.push_back(Triplet<double>(i, i + SIZE_J, 1));
			}
		}
	}
	/*SparseMatrix<double> matA(SIZE_I * SIZE_J, SIZE_I * SIZE_J);*/
	matA.setFromTriplets(tripletList.begin(), tripletList.end());
	std::cout << matA << std::endl;
	double x = matA.coeff(0, 0);
	std::cout << x << std::endl;
	cg.setMaxIterations(1);
	//cg.compute(matA);
}

// The projection step using conjugate gradient
void MyWorld::project_CG(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _d[][SIZE_J + 2]) {
	
	//////////////// Use the conjugate gradient to calculate pressure p ////////////////////////
	std::cout << "Divergence Vector" << std::endl;
	VectorXd divVect(SIZE_I * SIZE_J);
	for (int i = 1; i <= SIZE_I; i++) {
		for (int j = 1; j <= SIZE_J; j++) {
			divVect((i - 1) + (SIZE_I * (j - 1))) = _d[i][j];
			std::cout << divVect((i - 1) + (SIZE_I * (j - 1))) << std::endl;
		}
	}
	/*std::cout << "divergence vector" << std::endl;
	std::cout << divVect << std::endl;
	std::cout << " " << std::endl;*/
	std::cout << "Divergence Vector" << std::endl;
	std::cout << divVect(3) << std::endl;
	std::cout << divVect(4) << std::endl;
	std::cout << divVect(5) << std::endl;
	std::cout << " " << std::endl;

	VectorXd pressure(SIZE_I * SIZE_J);
	cg.setMaxIterations(2);
	pressure = cg.compute(matA).solve(divVect);

	std::cout << "Pressure Vector" << std::endl;
	std::cout << pressure(3) << std::endl;
	std::cout << pressure(4) << std::endl;
	std::cout << pressure(5) << std::endl;
	std::cout << " " << std::endl;
	
	/////////////// Find the gradient of pressure p /////////////////////////////////
	double h = (SIZE_I + SIZE_J) / 2.0;
	double pGradU[SIZE_I][SIZE_J];
	double pGradV[SIZE_I][SIZE_J];
	for (int i = 0; i < SIZE_I; i++) {
		for (int j = 0; j < SIZE_J; j++) {
			if (i == SIZE_I - 1) {
				pGradU[i][j] = 0.0;
			}
			else {
				pGradU[i][j] = (pressure((i + 1) * j) - pressure(i * j)) / h;
			}
			if (j == SIZE_I - 1) {
				pGradV[i][j] = 0.0;
			}
			else {
				pGradV[i][j] = (pressure(i * (j + 1)) - pressure(i * j)) / h;
			}
		}
	}

	/*std::cout << _u[40][40] << std::endl;
	std::cout << _v[40][40] << std::endl;*/

	/*std::cout << " " << std::endl;
	std::cout << " " << std::endl;*/

	//////////////// Update velocity from the gradient of pressure ///////////////////////////
	for (int i = 1; i <= SIZE_I; i++) {
		for (int j = 1; j <= SIZE_J; j++) {
			_u[i][j] -= pGradU[i - 1][j - 1];
			_v[i][j] -= pGradV[i - 1][j - 1];
		}
	}
	/*std::cout << _u[5][5] << std::endl;
	std::cout << _v[5][5] << std::endl;

	std::cout << " " << std::endl;
	std::cout << " " << std::endl;*/

}

double MyWorld::getTimeStep() {
    return mTimeStep;
}

void MyWorld::simulate() {
    //velocityStep_1D(mU, mV, mPreU, mPreV);
	velocityStep_2D(mU, mV, mPreU, mPreV);
	/*calcDivergence(mDivergenceArray);*/
    //densityStep(mDensity, mPreDensity);
	/*densityStep_1D(mDensity_R, mPreDensity_R);
	densityStep_1D(mDensity_G, mPreDensity_G);
	densityStep_1D(mDensity_B, mPreDensity_B);*/
	densityStep_2D(mDensity_R, mPreDensity_R);
	densityStep_2D(mDensity_G, mPreDensity_G);
	densityStep_2D(mDensity_B, mPreDensity_B);
    /*externalForces_1D();*/
	externalForces_2D();
}

void MyWorld::densityStep_1D(double *_x, double *_x0) {
    SWAPPING_1D(_x, _x0); // _x now points at mPreDensity
    diffuseDensity_1D(_x, _x0); // Diffusion on _x which points at mPreDensity
    SWAPPING_1D(_x, _x0); // _x now points at mDensity
    //advectDensity_1D(_x, _x0, mU, mV); // Advection on _x which points at mDensity
}

//void MyWorld::densityStep_2D(double **_x, double **_x0) {
//	SWAPPING_2D(_x, _x0); // _x now points at mPreDensity
//	diffuseDensity_2D(_x, _x0); // Diffusion on _x which points at mPreDensity
//	SWAPPING_2D(_x, _x0); // _x now points at mDensity
//	advectDensity_2D(_x, _x0, mU, mV); // Advection on _x which points at mDensity
//}


void MyWorld::densityStep_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2]) {
	//SWAPPING_2D(_x, _x0); // _x now points at mPreDensity
	diffuseDensity_2D(_x0, _x); // Diffusion on _x which points at mPreDensity
	//SWAPPING_2D(_x, _x0); // _x now points at mDensity
	advectDensity_2D(_x, _x0, mU, mV); // Advection on _x which points at mDensity
}


void MyWorld::velocityStep_1D(double *_u, double *_v, double *_u0, double *_v0) {
    SWAPPING_1D(_u, _u0); // _u now points at mPreU
    SWAPPING_1D(_v, _v0); // _v now points at mPreV
    diffuseVelocity_1D(_u, _v, _u0, _v0);
    project_1D(_u, _v, _u0, _v0);
    SWAPPING_1D(_u, _u0); // _u now points at mU
    SWAPPING_1D(_v, _v0); // _u now points at mV
    advectVelocity_1D(_u, _v, _u0, _v0);
    project_1D(_u, _v, _u0, _v0);
}

void MyWorld::velocityStep_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]) {
	if (frameCount < 4) {
		//SWAPPING_1D(_u, _u0); // _u now points at mPreU
		//SWAPPING_1D(_v, _v0); // _v now points at mPreV
		diffuseVelocity_2D(_u0, _v0, _u, _v);
		calcDivergence(mDivergenceArray);
		cout << "Frame Count: " << frameCount << "\n";
		frameCount++;
		cout << "(Before project) Average Divergence: " << averageDivergence << "\n";
		cout << "(Before project) Average Absolute Divergence: " << averageAbsDivergence << "\n";
		//project_2D(_u0, _v0, _u, _v);
		project_CG(_u, _v, mDivergenceArray);
		calcDivergence(mDivergenceArray);
		cout << "(After project) Average Divergence: " << averageDivergence << "\n";
		cout << "(After project) Average Absolute Divergence: " << averageAbsDivergence << "\n\n";
		//SWAPPING_1D(_u, _u0); // _u now points at mU
		//SWAPPING_1D(_v, _v0); // _u now points at mV
		advectVelocity_2D(_u, _v, _u0, _v0);
		//project_2D(_u, _v, _u0, _v0);
		project_CG(_u, _v, mDivergenceArray);
		calcDivergence(mDivergenceArray);
	}
	
}

void MyWorld::diffuseDensity_1D(double *_x, double *_x0) {
    double a = mTimeStep * mDiffusionCoef * mNumCells * mNumCells;
    linearSolve_1D(_x, _x0, a, 1 + 4 * a);
    setBoundary_1D(_x);
}

void MyWorld::diffuseDensity_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2]) {
	double a = mTimeStep * mDiffusionCoef * SIZE_I * SIZE_J;
	linearSolve_2D(_x, _x0, a, 1 + 4 * a);
	setBoundary_2D(_x);
}

void MyWorld::diffuseVelocity_1D(double *_u, double *_v, double *_u0, double *_v0) {
    double a = mTimeStep * mViscosityCoef * mNumCells * mNumCells;
    linearSolve_1D(_u, _u0, a, 1 + 4 * a);
    linearSolve_1D(_v, _v0, a, 1 + 4 * a);
    setVelocityBoundary_1D(_u, _v);
}

void MyWorld::diffuseVelocity_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]) {
	double a = mTimeStep * mViscosityCoef * SIZE_I * SIZE_J;
	linearSolve_2D(_u, _u0, a, 1 + 4 * a);
	linearSolve_2D(_v, _v0, a, 1 + 4 * a);
	setVelocityBoundary_2D(_u, _v);
}

void MyWorld::advectDensity_1D(double *_d, double *_d0, double *_u, double *_v) {
    double dt0 = mTimeStep * mNumCells;  // h / x 
    for (int i = 1; i <= mNumCells; i++) {
        for (int j = 1; j <= mNumCells; j++) {
          //          if (abs(_u[IX(i, j)]) > 0.001)
            
            double x = i- dt0 * _u[IX(i,j)];  // dt0 * _u[IX(i,j)] computes how many cells can a particle travel in one time step 
            double y = j - dt0 * _v[IX(i,j)];
            if (x < 0.5) 
                x = 0.5f; 
            if (x > mNumCells + 0.5) 
                x = mNumCells + 0.5;
            int i0 = (int)x;
            int i1 = i0 + 1;
            if (y < 0.5) 
                y = 0.5;
            if (y > mNumCells + 0.5)
                y = mNumCells + 0.5;
            int j0 = (int)y;
            int j1 = j0 + 1;
            double s1 = x - i0;
            double s0 = 1 - s1;
            double t1 = y - j0;
            double t0 = 1 - t1;
            _d[IX(i,j)] = s0 * (t0 * _d0[IX(i0, j0)] + t1 * _d0[IX(i0,j1)])+ s1 * (t0 * _d0[IX(i1, j0)] + t1 * _d0[IX(i1,j1)]);
		}
    }
    setBoundary_1D(_d);
}

//void MyWorld::advectDensity_2D(double _d[][SIZE + 2], double _d0[][SIZE + 2], double *_u, double *_v) {
//	double dt0 = mTimeStep * mNumCells;  // h / x 
//	for (int i = 1; i <= mNumCells; i++) {
//		for (int j = 1; j <= mNumCells; j++) {
//			//          if (abs(_u[IX(i, j)]) > 0.001)
//
//			double x = i - dt0 * _u[IX(i, j)];  // dt0 * _u[i][j] computes how many cells can a particle travel in one time step 
//			double y = j - dt0 * _v[IX(i, j)];
//			if (x < 0.5)
//				x = 0.5f;
//			if (x > mNumCells + 0.5)
//				x = mNumCells + 0.5;
//			int i0 = (int)x;
//			int i1 = i0 + 1;
//			if (y < 0.5)
//				y = 0.5;
//			if (y > mNumCells + 0.5)
//				y = mNumCells + 0.5;
//			int j0 = (int)y;
//			int j1 = j0 + 1;
//			double s1 = x - i0;
//			double s0 = 1 - s1;
//			double t1 = y - j0;
//			double t0 = 1 - t1;
//			_d[i][j] = s0 * (t0 * _d0[i0][j0] + t1 * _d0[i0][j1]) + s1 * (t0 * _d0[i1][j0] + t1 * _d0[i1][j1]);
//		}
//	}
//	setBoundary_2D(_d);
//}

void MyWorld::advectDensity_2D(double _d[][SIZE_J + 2], double _d0[][SIZE_J + 2], double _u[][SIZE_J + 2], double _v[][SIZE_J + 2]) {
	double dt0_x = mTimeStep * SIZE_I;  // h / x 
	double dt0_y = mTimeStep * SIZE_J;  // h / y 

	for (int i = 1; i <= SIZE_I; i++) {
		for (int j = 1; j <= SIZE_J; j++) {
			//          if (abs(_u[IX(i, j)]) > 0.001)

			double x = i - dt0_x * _u[i][j];  // dt0 * _u[i][j] computes how many cells can a particle travel in one time step 
			double y = j - dt0_y * _v[i][j];
			if (x < 0.5)
				x = 0.5f;
			if (x > SIZE_I + 0.5)
				x = SIZE_I + 0.5;
			int i0 = (int)x;
			int i1 = i0 + 1;
			if (y < 0.5)
				y = 0.5;
			if (y > SIZE_J + 0.5)
				y = SIZE_J + 0.5;
			int j0 = (int)y;
			int j1 = j0 + 1;
			double s1 = x - i0;
			double s0 = 1 - s1;
			double t1 = y - j0;
			double t0 = 1 - t1;
			_d[i][j] = s0 * (t0 * _d0[i0][j0] + t1 * _d0[i0][j1]) + s1 * (t0 * _d0[i1][j0] + t1 * _d0[i1][j1]);
		}
	}
	setBoundary_2D(_d);
}

void MyWorld::advectVelocity_1D(double *_u, double *_v, double *_u0, double *_v0) {
    // TODO: Add velocity advection code here
	double dt0 = mTimeStep * mNumCells;  // h / x 
	for (int i = 1; i <= mNumCells; i++) {
		for (int j = 1; j <= mNumCells; j++) {
			//          if (abs(_u[IX(i, j)]) > 0.001)

			double x = i - dt0 * _u0[IX(i, j)];  // dt0 * _u[IX(i,j)] computes how many cells can a particle travel in one time step 
			double y = j - dt0 * _v0[IX(i, j)];
			if (x < 0.5)
				x = 0.5f;
			if (x > mNumCells + 0.5)
				x = mNumCells + 0.5;
			int i0 = (int)x;
			int i1 = i0 + 1;
			if (y < 0.5)
				y = 0.5;
			if (y > mNumCells + 0.5)
				y = mNumCells + 0.5;
			int j0 = (int)y;
			int j1 = j0 + 1;
			double s1 = x - i0;
			double s0 = 1 - s1;
			double t1 = y - j0;
			double t0 = 1 - t1;
			_u[IX(i, j)] = s0 * (t0 * _u0[IX(i0, j0)] + t1 * _u0[IX(i0, j1)]) + s1 * (t0 * _u0[IX(i1, j0)] + t1 * _u0[IX(i1, j1)]);
			_v[IX(i, j)] = s0 * (t0 * _v0[IX(i0, j0)] + t1 * _v0[IX(i0, j1)]) + s1 * (t0 * _v0[IX(i1, j0)] + t1 * _v0[IX(i1, j1)]);
		}
	}
	setVelocityBoundary_1D(_u, _v);
}

void MyWorld::advectVelocity_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]) {
	// TODO: Add velocity advection code here
	double dt0_x = mTimeStep * SIZE_I;  // h / x 
	double dt0_y = mTimeStep * SIZE_J;  // h / y 

	for (int i = 1; i <= SIZE_I; i++) {
		for (int j = 1; j <= SIZE_J; j++) {
			//          if (abs(_u[IX(i, j)]) > 0.001)

			double x = i - dt0_x * _u0[i][j];  // dt0 * _u[i][j] computes how many cells can a particle travel in one time step 
			double y = j - dt0_y * _v0[i][j];
			if (x < 0.5)
				x = 0.5f;
			if (x > SIZE_I + 0.5)
				x = SIZE_I + 0.5;
			int i0 = (int)x;
			int i1 = i0 + 1;
			if (y < 0.5)
				y = 0.5;
			if (y > SIZE_J + 0.5)
				y = SIZE_J + 0.5;
			int j0 = (int)y;
			int j1 = j0 + 1;
			double s1 = x - i0;
			double s0 = 1 - s1;
			double t1 = y - j0;
			double t0 = 1 - t1;
			_u[i][j] = s0 * (t0 * _u0[i0][j0] + t1 * _u0[i0][j1]) + s1 * (t0 * _u0[i1][j0] + t1 * _u0[i1][j1]);
			_v[i][j] = s0 * (t0 * _v0[i0][j0] + t1 * _v0[i0][j1]) + s1 * (t0 * _v0[i1][j0] + t1 * _v0[i1][j1]);
		}
	}
	setVelocityBoundary_2D(_u, _v);
}

void MyWorld::project_1D(double *_u, double *_v, double *_u0, double *_v0) {
   // TODO: Add projection code here
	int i, j, k;
	double h;

	h = 1.0 / mNumCells;
	for (i = 1; i <= mNumCells; i++) {
		for (j = 1; j <= mNumCells; j++) {
			_v0[IX(i,j)] = -0.5 * h * (_u[IX(i + 1, j)] - _u[IX(i - 1, j)] + _v[IX(i, j + 1)] - _v[IX(i, j - 1)]);
			_u0[IX(i, j)] = 0;
		}
	}

	setBoundary_1D(_v0);
	setBoundary_1D(_u0);
	
	for (k = 0; k < 20; k++) {
		for (i = 1; i <= mNumCells; i++) {
			for (j = 1; j <= mNumCells; j++) {
				_u0[IX(i, j)] = (_v0[IX(i, j)] + _u0[IX(i - 1, j)] + _u0[IX(i + 1, j)] + _u0[IX(i, j - 1)] + _u0[IX(i, j + 1)]) / 4.0;
			}
		}
		setBoundary_1D(_u0);
	}
	//linearSolve(_u0, _v0, 1, 4);

	for (i = 1; i <= mNumCells; i++) {
		for (j = 1; j <= mNumCells; j++) {
			_u[IX(i, j)] -= 0.5 * (_u0[IX(i + 1, j)] - _u0[IX(i - 1, j)]) / h;
			_v[IX(i, j)] -= 0.5 * (_u0[IX(i, j + 1)] - _u0[IX(i, j - 1)]) / h;
		}
	}

	setVelocityBoundary_1D(_u, _v);
}

void MyWorld::project_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2], double _u0[][SIZE_J + 2], double _v0[][SIZE_J + 2]) {
	// TODO: Add projection code here
	int i, j, k;
	double hx;
	double hy;
	double h;

	hx = 1.0 / SIZE_I;
	hy = 1.0 / SIZE_J;
	h = 1.0 / ((SIZE_I + SIZE_J) / 2.0);
	for (i = 1; i <= SIZE_I; i++) {
		for (j = 1; j <= SIZE_J; j++) {
			_v0[i][j] = -0.5 * h * (_u[i + 1][j] - _u[i - 1][j] + _v[i][j + 1] - _v[i][j - 1]);
			_u0[i][j] = 0;
		}
	}

	setBoundary_2D(_v0);
	setBoundary_2D(_u0);

	for (k = 0; k < 20; k++) {
		for (i = 1; i <= SIZE_I; i++) {
			for (j = 1; j <= SIZE_J; j++) {
				_u0[i][j] = (_v0[i][j] + _u0[i - 1][j] + _u0[i + 1][j] + _u0[i][j - 1] + _u0[i][j + 1]) / 4.0;
			}
		}
		setBoundary_2D(_u0);
	}
	//linearSolve(_u0, _v0, 1, 4);

	for (i = 1; i <= SIZE_I; i++) {
		for (j = 1; j <= SIZE_J; j++) {
			_u[i][j] -= 0.5 * (_u0[i + 1][j] - _u0[i - 1][j]) / h;
			_v[i][j] -= 0.5 * (_u0[i][j + 1] - _u0[i][j - 1]) / h;
		}
	}

	setVelocityBoundary_2D(_u, _v);
}

//void MyWorld::externalForces_1D() {
//    int size = (mNumCells + 2) * (mNumCells + 2);
//    for (int i = 0; i< size; i++) {
//        //mPreDensity[i] = 0;
//		mPreDensity_R[i] = 0;
//		mPreDensity_G[i] = 0;
//		mPreDensity_B[i] = 0;
//        mPreU[i] = 0;
//        mPreV[i] = 0;
//    }
//}

void MyWorld::externalForces_2D() {
	//int index = 0;
	for (int i = 0; i < sizeI; i++) {
		for (int j = 0; j < sizeJ; j++) {
			mPreDensity_R[i][j] = 0;
			mPreDensity_G[i][j] = 0;
			mPreDensity_B[i][j] = 0;

			/*mPreU[index] = 0;
			mPreV[index] = 0;*/
			//index++;
			mPreU[i][j] = 0;
			mPreV[i][j] = 0;
		}
	}
}

void MyWorld::calcDivergence(double _d[][SIZE_J + 2]) {
	double divergenceSum = 0.0;
	double absDivergenceSum = 0.0;
	double divergence = 0.0;

	//double scale = 1.0 / ((SIZE_I + SIZE_J) / 2.0);
	double scale = 1.0;

	for (int i = 1; i <= SIZE_I; i++) {
		for (int j = 1; j <= SIZE_J; j++) {
			if (j == 1) {
				if (i == 1) {
					divergence = mU[2][1] + mV[1][2];
				}
				else if (i == SIZE_I) {
					divergence = (-1.0 * mU[SIZE_I - 1][1]) + mV[SIZE_I][2];
				}
				else {
					divergence = (mU[i + 1][1] - mU[i - 1][1]) + mV[i][2];
				}
			} else if (j == SIZE_J) {
				if (i == 1) {
					divergence = mU[2][SIZE_J] - mV[1][SIZE_J - 1];
				}
				else if (i == SIZE_I) {
					divergence = (-1.0 * mU[SIZE_I - 1][SIZE_J]) - mV[SIZE_I][SIZE_J - 1];
				}
				else {
					divergence = (mU[i + 1][SIZE_J] - mU[i - 1][SIZE_J]) - mV[i][SIZE_J - 1];
				}
			} else {
				if (i == 1) {
					divergence = mU[1][j] + (mV[1][j + 1] - mV[1][j - 1]);
				}
				else if (i == SIZE_I) {
					divergence = (-1.0 * mU[SIZE_I - 1][j]) + (mV[SIZE_I][j + 1] - mV[SIZE_I][j - 1]);
				}
				else {
					divergence = (mU[i + 1][j] - mU[i - 1][j]) + (mV[i][j + 1] - mV[i][j - 1]);
				}
			}
			//divergence = -1.0 * scale * divergence;
			_d[i][j] = divergence;
			divergenceSum += divergence;
			absDivergenceSum += abs(divergence);
		}
	}
	averageDivergence = divergenceSum / (SIZE_I * SIZE_J);
	averageAbsDivergence = absDivergenceSum / (SIZE_I * SIZE_J);

	/*std::cout << "Average Absolute Divergence:" << std::endl;
	std::cout << averageAbsDivergence << std::endl;
	std::cout << "" << std::endl;*/
}

void MyWorld::linearSolve_1D(double *_x, double *_x0, double _a, double _c) {
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= mNumCells; i++) {
            for (int j = 1; j <= mNumCells; j++) {
                _x[IX(i, j)] = (_x0[IX(i, j)] + _a * (_x[IX(i-1, j)] + _x[IX(i+1, j)] + _x[IX(i, j-1)] + _x[IX(i, j+1)])) / _c;
            }
        }
    }
}

void MyWorld::linearSolve_2D(double _x[][SIZE_J + 2], double _x0[][SIZE_J + 2], double _a, double _c) {
	for (int k = 0; k < 20; k++) {
		for (int i = 1; i <= SIZE_I; i++) {
			for (int j = 1; j <= SIZE_J; j++) {
				_x[i][j] = (_x0[i][j] + _a * (_x[i - 1][j] + _x[i + 1][j] + _x[i][j - 1] + _x[i][j + 1])) / _c;
			}
		}
	}
}

void MyWorld::setBoundary_1D(double *_x) {
    for (int i = 1; i <= mNumCells; i++) {
        _x[IX(0 ,i)] = _x[IX(1,i)];
        _x[IX(mNumCells+1, i)] = _x[IX(mNumCells, i)];
        _x[IX(i, 0)] = _x[IX(i, 1)];
        _x[IX(i, mNumCells+1)] = _x[IX(i, mNumCells)];
 
    }
    _x[IX(0, 0)] = 0.5 * (_x[IX(1, 0)] + _x[IX(0, 1)]);
    _x[IX(0, mNumCells+1)] = 0.5 * (_x[IX(1, mNumCells+1)] + _x[IX(0, mNumCells)]);
    _x[IX(mNumCells+1, 0)] = 0.5 * (_x[IX(mNumCells, 0)] + _x[IX(mNumCells+1, 1)]);
    _x[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_x[IX(mNumCells, mNumCells+1)] + _x[IX(mNumCells+1, mNumCells)]);
}

void MyWorld::setBoundary_2D(double _x[][SIZE_J + 2]) {
	for (int i = 1; i <= SIZE_I; i++) {
		_x[0][i] = _x[1][i];
		_x[SIZE_I + 1][i] = _x[SIZE_I][i];
		_x[i][0] = _x[i][1];
		_x[i][SIZE_J + 1] = _x[i][SIZE_J];

	}
	_x[0][0] = 0.5 * (_x[1][0] + _x[0][1]);
	_x[0][SIZE_J + 1] = 0.5 * (_x[1][SIZE_J + 1] + _x[0][SIZE_J]);
	_x[SIZE_I + 1][0] = 0.5 * (_x[SIZE_I][0] + _x[SIZE_I + 1][1]);
	_x[SIZE_I + 1][SIZE_J + 1] = 0.5 * (_x[SIZE_I][SIZE_J + 1] + _x[SIZE_I + 1][SIZE_J]);
}

void MyWorld::setVelocityBoundary_1D(double *_u, double *_v) {
    for (int i = 1; i <= mNumCells; i++) {
        _u[IX(0 ,i)] = -_u[IX(1,i)];
        _u[IX(mNumCells+1, i)] = -_u[IX(mNumCells, i)];
        _u[IX(i, 0)] = _u[IX(i, 1)];
        _u[IX(i, mNumCells+1)] = _u[IX(i, mNumCells)];

        _v[IX(0 ,i)] = _v[IX(1,i)];
        _v[IX(mNumCells+1, i)] = _v[IX(mNumCells, i)];
        _v[IX(i, 0)] = -_v[IX(i, 1)];
        _v[IX(i, mNumCells+1)] = -_v[IX(i, mNumCells)];
    }
    _u[IX(0, 0)] = 0.5 * (_u[IX(1, 0)] + _u[IX(0, 1)]);
    _u[IX(0, mNumCells+1)] = 0.5 * (_u[IX(1, mNumCells+1)] + _u[IX(0, mNumCells)]);
    _u[IX(mNumCells+1, 0)] = 0.5 * (_u[IX(mNumCells, 0)] + _u[IX(mNumCells+1, 1)]);
    _u[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_u[IX(mNumCells, mNumCells+1)] + _u[IX(mNumCells+1, mNumCells)]);
    _v[IX(0, 0)] = 0.5 * (_v[IX(1, 0)] + _v[IX(0, 1)]);
    _v[IX(0, mNumCells+1)] = 0.5 * (_v[IX(1, mNumCells+1)] + _v[IX(0, mNumCells)]);
    _v[IX(mNumCells+1, 0)] = 0.5 * (_v[IX(mNumCells, 0)] + _v[IX(mNumCells+1, 1)]);
    _v[IX(mNumCells+1, mNumCells+1)] = 0.5 * (_v[IX(mNumCells, mNumCells+1)] + _v[IX(mNumCells+1, mNumCells)]);

}

void MyWorld::setVelocityBoundary_2D(double _u[][SIZE_J + 2], double _v[][SIZE_J + 2]) {
	for (int i = 1; i <= SIZE_I; i++) {
		_u[0][i] = -_u[1][i];
		_u[SIZE_I + 1][i] = -_u[SIZE_I][i];
		_u[i][0] = _u[i][1];
		_u[i][SIZE_J + 1] = _u[i][SIZE_J];

		_v[0][i] = _v[1][i];
		_v[SIZE_I + 1][i] = _v[SIZE_I][i];
		_v[i][0] = -_v[i][1];
		_v[i][SIZE_J + 1] = -_v[i][SIZE_J];
	}
	_u[0][0] = 0.5 * (_u[1][0] + _u[0][1]);
	_u[0][SIZE_J + 1] = 0.5 * (_u[1][SIZE_J + 1] + _u[0][SIZE_J]);
	_u[SIZE_I + 1][0] = 0.5 * (_u[SIZE_I][0] + _u[SIZE_I + 1][1]);
	_u[SIZE_I + 1][SIZE_J + 1] = 0.5 * (_u[SIZE_I][SIZE_J + 1] + _u[SIZE_I + 1][SIZE_J]);

	_v[0][0] = 0.5 * (_v[1][0] + _v[0][1]);
	_v[0][SIZE_J + 1] = 0.5 * (_v[1][SIZE_J + 1] + _v[0][SIZE_J]);
	_v[SIZE_I + 1][0] = 0.5 * (_v[SIZE_I][0] + _v[SIZE_I + 1][1]);
	_v[SIZE_I + 1][SIZE_J + 1] = 0.5 * (_v[SIZE_I][SIZE_J + 1] + _v[SIZE_I + 1][SIZE_J]);

}

//void MyWorld::swap_2D(double _x[][], double _x0[][64]) {
//	
//}
