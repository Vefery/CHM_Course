#pragma once
#include "stdafx.h"

#define REALOUTD "%.2f\t"

class SLAE {
public:
	int n, maxiter = 10000, nProfile = 0;
	double  eps = 1e-13;
	double* al, * au, * di;
	double* alLU, * auLU, * diLU;
	double* x, * x0, * b, * xtrue, * dP;
	double* r, * z, * tmp1;
	int* ia, * ja;

	void Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf);
	void Input(int n, int maxiter, double eps, int* ig, int* jg, double* gglA, double* gguA, double* diA, double* b);

	void MatrixVectorMultiplication(double* vectorMult, double* vectorOut) const;
	void TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut) const;
	double CalculateRelativeDiscrepancyWithR(double norm) const;
	double CalculateRelativeDiscrepancy(double norm) const;

	void MethodOfConjugateGradientsForSymMatrix();
	void MethodOfConjugateGradientsForNonSymMatrix();

	void MethodOfConjugateGradientsForSymMatrixWithDiagP();
	void MethodOfConjugateGradientsForNonSymMatrixWithDiagP();

	void MethodOfConjugateGradientsForSymMatrixWithLuP();
	void MethodOfConjugateGradientsForNonSymMatrixWithLuP();
	void MethodOfConjugateGradientsForNonSymMatrixWithLuAsterP();
	void MethodOfConjugateGradientsForNonSymMatrixWithLuSqP();

	void VectorConditionalityForSymMatrixDiagP(double* vectorIn, double* vectorOut) const;
	void VectorConditionalityForNonSymMatrixDiagP(double* vectorIn, double* vectorOut) const;

	void CalculateLU();
	void CalculateLUaster();
	void CalculateLUsq();

	void GenerateHilbertMatrix(int size);

	void SolveForwardLU(double* lowerTringMat, double* rightVector, double* vectorX) const;
	void SolveBackwardLU(double* upperTringMat, double* rightVector, double* vectorX) const;
	void SolveForwardLU(double* lowerTringMat, double* diag, double* rightVector, double* vectorX) const;
	void SolveBackwardLU(double* upperTringMat, double* diag, double* rightVector, double* vectorX) const;

	void MatrixUVectorMultiplicationLU(double* U, double* vectorMult, double* vectorOut) const;
	void MatrixUVectorMultiplicationLU(double* U, double *diag, double* vectorMult, double* vectorOut) const;
	void VectorUVectorMultiplication(double* x) const;
	void CalculateZ_LU(double* vectorOut) const;
	void CalculateZ_LUaster(double* vectorOut) const;
	void CalculateZ_LUsq(double* vectorOut) const;

	void CalculateXkRk(double ak);
	void CalculateZk(double bk);
	void Calculate_dP();

	void CalculateRelativeDiscrepancy(double* vectorMult, double* vectorOut) const;
	void VectorSubtract(double* first, double* second, double* result) const;
	double VectorScalarProduction(double* vector) const;
	double VectorScalarProduction(double* vector1, double* vector2) const;
	double VectorNorm(double* vector) const;
	void VectorCopy(double* first, double* second) const;

	void OutputDense();
	void OutputLUDense();


	void VectorOutputSolution(FILE* out);

protected:
	void AllocateMemory();
	void ClearMemory();
};