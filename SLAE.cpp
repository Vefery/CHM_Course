#include "stdafx.h"
#include "slae.h"

void SLAE::Input(FILE* paramf, FILE* iaf, FILE* jaf, FILE* alf, FILE* auf, FILE* dif, FILE* bf) {
	fscanf_s(paramf, "%d", &n);
	fscanf_s(paramf, "%d", &maxiter);
	fscanf_s(paramf, "%lf", &eps);

	ia = new int[n + 1];
	for (int i = 0; i <= n; i++)
		fscanf_s(iaf, "%d", &ia[i]);
	nProfile = ia[n] - ia[0];
	AllocateMemory();
	for (int i = 0; i < nProfile; i++)
		fscanf_s(jaf, "%d", &ja[i]);
	if (ia[0]) {
		for (int i = 0; i <= n; i++)
			ia[i]--;
		for (int i = 0; i < nProfile; i++)
			ja[i]--;
	}

	for (int i = 0; i < nProfile; i++)
		fscanf_s(alf, "%lf", &al[i]);

	for (int i = 0; i < nProfile; i++)
		fscanf_s(auf, "%lf", &au[i]);

	for (int i = 0; i < n; i++)
		fscanf_s(dif, "%lf", &di[i]);

	for (int i = 0; i < n; i++)
		fscanf_s(bf, "%lf", &b[i]);

	//for (int i = 0; i < n; i++)
	//	x0[i] = 0;
	for (int i = 0; i < n; i++)
	{
		xtrue[i] = i + 1;
	}
}

void SLAE::Input(int n, int maxiter, double eps, int* ig, int* jg, double* ggl, double* ggu, double* di, double* b)
{
	this->n = n;
	this->maxiter = maxiter;
	this->eps = eps;

	ia = new int[n + 1];
	std::copy(ig, ig + n + 1, ia);
	nProfile = ia[n] - ia[0];
	AllocateMemory();
	std::copy(jg, jg + nProfile, ja);

	if (ia[0]) {
		for (int i = 0; i <= n; i++)
			ia[i]--;
		for (int i = 0; i < nProfile; i++)
			ja[i]--;
	}

	std::copy(ggl, ggl + nProfile, al);

	std::copy(ggu, ggu + nProfile, au);

	std::copy(di, di + n, this->di);
	
	std::copy(b, b + n, this->b);

	for (int i = 0; i < n; i++)
	{
		xtrue[i] = i + 1;
	}
}

void SLAE::OutputDense()
{
	int flagfound = 0;
	for (int i = 0; i < n; i++)
	{
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++) //��� ��� 100 ���� ���������.
			{
				for (int p = lastj; p < ja[j]; p++) //��� ��� 100 ���� ���������.
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, al[j]);
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++) //??
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, di[i]);

		for (int j = i + 1; j < n; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, au[k]);
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}
	printf("\n");
}

void SLAE::OutputLUDense()
{
	int flagfound = 0;
	for (int i = 0; i < n; i++)
	{
		int k = ia[i + 1] - ia[i];
		if (k == 0)
		{
			for (int j = 0; j < i; j++)
			{
				printf(REALOUTD, 0.0);
			}
		}
		else
		{
			int lastj = 0;
			for (int j = ia[i]; j < ia[i + 1]; j++) //��� ��� 100 ���� ���������.
			{
				for (int p = lastj; p < ja[j]; p++) //��� ��� 100 ���� ���������.
				{
					printf(REALOUTD, 0.0);
				}
				printf(REALOUTD, alLU[j]);
				lastj = ja[j] + 1;
			}
			for (int j = lastj; j < i; j++) //??
			{
				printf(REALOUTD, 0.0);
			}
		}

		printf(REALOUTD, diLU[i]);

		for (int j = i + 1; j < n; j++)
		{
			k = ia[j + 1] - ia[j];
			if (k == 0) {
				printf(REALOUTD, 0.0);
			}
			else
			{
				flagfound = 0;
				for (k = ia[j]; k < ia[j + 1]; k++)
				{

					if (ja[k] == i)
					{
						printf(REALOUTD, auLU[k]);
						flagfound = 1;
						break;
					}
				}
				if (flagfound == 0)
					printf(REALOUTD, 0.0);
			}
		}
		printf("\n");
	}
	printf("\n");
}

void SLAE::MatrixVectorMultiplication(double* vectorMult, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += al[k] * vectorMult[j];
			vectorOut[j] += au[k] * vectorMult[i];
		}
	}
}

void SLAE::TransposedMatrixVectorMultiplication(double* vectorMult, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += au[k] * vectorMult[j];
			vectorOut[j] += al[k] * vectorMult[i];
		}
	}
}

//*********************** ������������ ������� ***********************************************

void SLAE::MethodOfConjugateGradientsForSymMatrix()
{
	VectorCopy(x0, x);
	CalculateRelativeDiscrepancy(x0, tmp1);
	VectorCopy(tmp1, r);
	VectorCopy(r, z);

	double normB = VectorNorm(b);
	double ak = 0, bk = 0;
	double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		MatrixVectorMultiplication(z, tmp1); // A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		CalculateXkRk(ak);

		double r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;

		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
}

void SLAE::MethodOfConjugateGradientsForSymMatrixWithDiagP()
{
	VectorCopy(x0, x);
	CalculateRelativeDiscrepancy(x0, tmp1);
	VectorCopy(tmp1, r);
	VectorConditionalityForSymMatrixDiagP(tmp1, z);



	VectorConditionalityForSymMatrixDiagP(r, tmp1);
	double r_rPrev = VectorScalarProduction(tmp1, r); //( M^(-1)*r(k-1) , r(k-1) )

	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);
	double ak = 0, bk = 0;
	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		//VectorConditionalityForSymMatrixDiagP(r, tmp1);
		//double r_rPrev = VectorScalarProduction(tmp1, r); //( M^(-1)*r(k-1) , r(k-1) )
		MatrixVectorMultiplication(z, tmp1); //A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp1, z); // ( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		CalculateXkRk(ak);

		VectorConditionalityForSymMatrixDiagP(r, tmp1); // M^(-1)*rk
		double r_rCur = VectorScalarProduction(tmp1, r); //( M^(-1)*rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;
	}
}

void SLAE::MethodOfConjugateGradientsForSymMatrixWithLuP()
{
	VectorCopy(x0, x);
	CalculateRelativeDiscrepancy(x0, r);
	// M = LU
	// M^(-1) = U^(-1)*L^(-1)
	SolveForwardLU(alLU, diLU, r, tmp1);
	SolveBackwardLU(auLU, tmp1, z);


	double normB = VectorNorm(b);
	double ak = 0, bk = 0;

	//ak
	SolveForwardLU(alLU, diLU, r, tmp1);
	SolveBackwardLU(auLU, tmp1, tmp1);
	double Mr_rPrev = VectorScalarProduction(tmp1, r); // ( M^(-1)*r(k-1) , r(k-1) )
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		MatrixVectorMultiplication(z, tmp1); // A*z(k-1)
		double Az_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = Mr_rPrev / Az_zPrev;

		CalculateXkRk(ak);

		SolveForwardLU(alLU, diLU, r, tmp1);
		SolveBackwardLU(auLU, tmp1, tmp1);
		double Mr_rCur = VectorScalarProduction(tmp1, r); // ( rk , rk )

		bk = Mr_rCur / Mr_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		Mr_rPrev = Mr_rCur;
		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
}

//*********************************************************************

void SLAE::MethodOfConjugateGradientsForNonSymMatrix()
{
	VectorCopy(x0, x);
	CalculateRelativeDiscrepancy(x0, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, r);
	VectorCopy(r, z);

	double normB = VectorNorm(b);
	double ak = 0, bk = 0;
	double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double r_rCur, Az_zPrev;
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		MatrixVectorMultiplication(z, x0); // A*z(k-1)
		TransposedMatrixVectorMultiplication(x0, tmp1);
		Az_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		CalculateXkRk(ak);

		r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);
		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);
		
		r_rPrev = r_rCur;
		
		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));
}

void SLAE::MethodOfConjugateGradientsForNonSymMatrixWithDiagP()
{
	Calculate_dP();

	VectorConditionalityForNonSymMatrixDiagP(x0, x);

	VectorCopy(x0, x);

	VectorUVectorMultiplication(x);

	CalculateRelativeDiscrepancy(x0, tmp1);

	VectorConditionalityForNonSymMatrixDiagP(tmp1, tmp1);
	VectorConditionalityForNonSymMatrixDiagP(tmp1, tmp1);

	TransposedMatrixVectorMultiplication(tmp1, x0);

	VectorConditionalityForNonSymMatrixDiagP(x0, x0);
	VectorCopy(x0, r);
	VectorCopy(x0, z);

	double normB = VectorNorm(b);
	double ak = 0, bk = 0;
	double r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double r_rCur, Az_zPrev;
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		VectorConditionalityForNonSymMatrixDiagP(z, tmp1);

		MatrixVectorMultiplication(tmp1, x0); // A*z(k-1)

		VectorConditionalityForNonSymMatrixDiagP(x0, x0); 
		VectorConditionalityForNonSymMatrixDiagP(x0, x0); 

		TransposedMatrixVectorMultiplication(x0, tmp1); // A^T*A*z(k - 1)

		VectorConditionalityForNonSymMatrixDiagP(tmp1, tmp1);

	    Az_zPrev = VectorScalarProduction(tmp1, z); // ( A^T*A*z(k-1) , z(k-1) )
		ak = r_rPrev / Az_zPrev;

		CalculateXkRk(ak);

	    r_rCur = VectorScalarProduction(r); // ( rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;

		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}

	VectorConditionalityForNonSymMatrixDiagP(x, x);

	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));
}

void SLAE::MethodOfConjugateGradientsForNonSymMatrixWithLuP()
{
	MatrixUVectorMultiplicationLU(auLU, x0, x);

	CalculateRelativeDiscrepancy(x0, tmp1);
	SolveForwardLU(alLU, diLU, tmp1, tmp1);
	SolveBackwardLU(alLU, diLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, x0, tmp1);
	VectorCopy(tmp1, r);
	VectorCopy(tmp1, z);

	double r_rPrev, r_rCur, Newz_zPrev, ak, bk;
	r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		VectorCopy(z, tmp1);
		CalculateZ_LU(tmp1); // U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
		Newz_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Newz_zPrev;

		CalculateXkRk(ak);

		r_rCur = VectorScalarProduction(r); // ( M^(-1)*rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;
		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
	SolveBackwardLU(auLU, x, x);
	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));
}

void SLAE::MethodOfConjugateGradientsForNonSymMatrixWithLuAsterP()
{
	MatrixUVectorMultiplicationLU(auLU, diLU, x0, x);

	CalculateRelativeDiscrepancy(x0, tmp1);
	SolveForwardLU(alLU, tmp1, tmp1);
	SolveBackwardLU(alLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, diLU, x0, tmp1);
	VectorCopy(tmp1, r);
	VectorCopy(tmp1, z);

	double r_rPrev, r_rCur, Newz_zPrev, ak, bk;
	r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		VectorCopy(z, tmp1);
		CalculateZ_LUaster(tmp1); // U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
		Newz_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Newz_zPrev;

		
		CalculateXkRk(ak);

		r_rCur = VectorScalarProduction(r); // ( M^(-1)*rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);

		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;

		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
	SolveBackwardLU(auLU, diLU, x, x);
	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));
}

void SLAE::MethodOfConjugateGradientsForNonSymMatrixWithLuSqP()
{
	MatrixUVectorMultiplicationLU(auLU, x0, x);

	CalculateRelativeDiscrepancy(x0, tmp1);
	SolveForwardLU(alLU, diLU, tmp1, tmp1);
	SolveBackwardLU(alLU, diLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, diLU, x0, tmp1);
	VectorCopy(tmp1, r);
	VectorCopy(tmp1, z);

	double r_rPrev, r_rCur, Newz_zPrev, ak, bk;
	r_rPrev = VectorScalarProduction(r); // ( r(k-1) , r(k-1) )
	double normB = VectorNorm(b);
	double RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

	for (int curIt = 0; curIt < maxiter and RelDiscrepancy > eps; curIt++)
	{
		VectorCopy(z, tmp1);
		CalculateZ_LUsq(tmp1); // U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
		Newz_zPrev = VectorScalarProduction(tmp1, z); //( A*z(k-1) , z(k-1) )
		ak = r_rPrev / Newz_zPrev;

		CalculateXkRk(ak);

		r_rCur = VectorScalarProduction(r); // ( M^(-1)*rk , rk )
		bk = r_rCur / r_rPrev;

		CalculateZk(bk);
		RelDiscrepancy = CalculateRelativeDiscrepancyWithR(normB);

		r_rPrev = r_rCur;

		printf("Iteration: %d, RelDiscrepancy of r: %.15lf\n", curIt + 1, RelDiscrepancy);
	}
	SolveBackwardLU(auLU, diLU, x, x);
	printf("%.15lf\n", CalculateRelativeDiscrepancy(normB));
}
//**********************************************************************


void SLAE::VectorUVectorMultiplication(double *x)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = x[i] * dP[i];
	}
}

void SLAE::CalculateXkRk(double ak)
{
	for (int i = 0; i < n; i++)
	{
		x[i] = x[i] + ak * z[i]; // xk = x(k-1) + ak*z(k-1)
		r[i] = r[i] - ak * tmp1[i]; // rk = r(k-1) - ak*U^(-T)*A^T*L^(-T)*L^(-1)*A*U^(-1)*z(k-1)
	}
}

void SLAE::CalculateZk(double bk)
{
	for (int i = 0; i < n; i++) // zk =  M^(-1)*rk + bk*z(k-1)
		z[i] = r[i] + bk * z[i];
}

void SLAE::CalculateLU()
{
	for (int i = 0; i < n; i++) //i = 7
	{
		int i0 = ia[i]; //15
		int i1 = ia[i + 1]; //19
		double sumD = 0;

		for (int k = i0; k < i1; k++)
		{
			double sumL = 0, sumU = 0;
			int j = ja[k];

			int j0 = ia[j];
			int j1 = ia[j + 1];

			int ik = i0;
			int kj = j0;

			for (; ik < k and kj < j1; )
			{
				if (ja[ik] < ja[kj])
					ik++;
				if (ja[ik] > ja[kj])
					kj++;
				if (ja[ik] == ja[kj]) {
					sumL += alLU[ik] * auLU[kj];
					sumU += alLU[kj] * auLU[ik];
					ik++; kj++;
				}

			}
			alLU[k] = al[k] - sumL;
			auLU[k] = (au[k] - sumU) / diLU[j];
			sumD += alLU[k] * auLU[k];
		}
		diLU[i] = di[i] - sumD;
	}
}

void SLAE::CalculateLUaster()
{
	for (int i = 0; i < n; i++) //i = 7
	{
		int i0 = ia[i]; //15
		int i1 = ia[i + 1]; //19
		double sumD = 0;

		for (int k = i0; k < i1; k++)
		{
			double sumL = 0, sumU = 0;
			int j = ja[k];

			int j0 = ia[j];
			int j1 = ia[j + 1];

			int ik = i0;
			int kj = j0;

			for (; ik < k and kj < j1; )
			{
				if (ja[ik] < ja[kj])
					ik++;
				if (ja[ik] > ja[kj])
					kj++;
				if (ja[ik] == ja[kj]) {
					sumL += alLU[ik] * auLU[kj];
					sumU += alLU[kj] * auLU[ik];
					ik++; kj++;
				}

			}
			alLU[k] = (al[k] - sumL) / diLU[j];
			auLU[k] = (au[k] - sumU);
			sumD += alLU[k] * auLU[k];
		}
		diLU[i] = di[i] - sumD;
	}
}

void SLAE::CalculateLUsq()
{
	for (int i = 0; i < n; i++) //i = 7
	{
		int i0 = ia[i]; //15
		int i1 = ia[i + 1]; //19
		double sumD = 0;

		for (int k = i0; k < i1; k++)
		{
			double sumL = 0, sumU = 0;
			int j = ja[k];

			int j0 = ia[j];
			int j1 = ia[j + 1];

			int ik = i0;
			int kj = j0;

			for (; ik < k and kj < j1; )
			{
				if (ja[ik] < ja[kj])
					ik++;
				if (ja[ik] > ja[kj])
					kj++;
				if (ja[ik] == ja[kj]) {
					sumL += alLU[ik] * auLU[kj];
					sumU += alLU[kj] * auLU[ik];
					ik++; kj++;
				}

			}
			alLU[k] = (al[k] - sumL) / diLU[j];
			auLU[k] = (au[k] - sumU) / diLU[j];
			sumD += alLU[k] * auLU[k];
		}
		diLU[i] = sqrt(di[i] - sumD);
	}
}

//**********************************************************************

void SLAE::GenerateHilbertMatrix(int size)
{
	n = size;
	for (int i = 0; i < size; i++)
	{
		nProfile += i;
	}
	ia = new int[n + 1];

	AllocateMemory();
	ia[0] = 0;
	for (int i = 1, k = 0; i <= n; i++)
	{
		ia[i] = ia[i - 1] + (i - 1);


		di[i - 1] = (double)1 / (2 * i - 1);
		for (int j = 1; j < i; j++, k++)
		{
			al[k] = (double)1 / (i + j - 1);
			au[k] = (double)1 / (i + j - 1);
			ja[k] = j - 1;
		}
	}

	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int xk = 1; xk <= n; xk++)
		{
			sum += (double)1 / (i + xk) * xk;
		}
		b[i] = sum;
	}
}


//**********************************************************************


void SLAE::CalculateRelativeDiscrepancy(double * vectorMult, double * vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = di[i] * vectorMult[i];
		vectorOut[i] = b[i] - vectorOut[i];
		for (int k = ia[i]; k < ia[i + 1]; k++)
		{
			int j = ja[k];
			vectorOut[i] += al[k] * vectorMult[j];
			vectorOut[j] += au[k] * vectorMult[i];
		}
	}
}

void SLAE::VectorConditionalityForSymMatrixDiagP(double* vectorIn, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorIn[i] / di[i];
	}
}

void SLAE::Calculate_dP()
{
	for (int i = 0; i < n; i++)
	{
		dP[i] = sqrt(di[i]);
	}
}

void SLAE::VectorConditionalityForNonSymMatrixDiagP(double* vectorIn, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorIn[i] / dP[i];
	}
}

void SLAE::VectorSubtract(double* first, double* second, double* result)
{
	for (int i = 0; i < n; i++)
	{
		result[i] = first[i] - second[i];
	}
}


void SLAE::VectorCopy(double* from, double* to)
{
	for (int i = 0; i < n; i++)
	{
		to[i] = from[i];
	}
}

double SLAE::VectorScalarProduction(double* vector1, double* vector2)
{
	double prod = 0;
	for (int i = 0; i < n; i++)
	{
		prod += vector1[i] * vector2[i];
	}
	return prod;
}

double SLAE::VectorScalarProduction(double* vector)
{
	double prod = 0;
	for (int i = 0; i < n; i++)
	{
		prod += vector[i] * vector[i];
	}
	return prod;
}

double SLAE::VectorNorm(double* vector)
{
	double norm = 0;
	for (int i = 0; i < n; i++)
	{
		norm += vector[i] * vector[i];
	}
	return sqrt(norm);
}

void SLAE::VectorOutputSolution(FILE* out)
{
	for (int i = 0; i < n; i++)
	{
		fprintf(out, "%.15lf\n", x[i]);
	}
	printf("\n");
	ClearMemory();
}


void SLAE::SolveForwardLU(double* lowerTringMat,double *diag, double* rightVector, double* vectorX) {
	for (int i0, i1, i = 0; i < n; i++)
	{
		double sum = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
		//int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++)
		{
			int j = ja[k];
			sum += lowerTringMat[k] * vectorX[j];
		}
		vectorX[i] = (rightVector[i] - sum) / diag[i];
	}
}

void SLAE::SolveForwardLU(double* lowerTringMat, double* rightVector, double* vectorX) {
	for (int i0, i1, i = 0; i < n; i++)
	{
		double sum = 0;
		i0 = ia[i];
		i1 = ia[i + 1];
		//int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++)
		{
			int j = ja[k];
			sum += lowerTringMat[k] * vectorX[j];
		}
		vectorX[i] = (rightVector[i] - sum);
	}
}


void SLAE::SolveBackwardLU(double* upperTringMat, double* diag, double* rightVector, double* vectorX) {
	for (int i0, i1, i = n - 1; i >= 0; i--)
	{
		i0 = ia[i];
		i1 = ia[i + 1];
		vectorX[i] = rightVector[i] / diag[i];
		for (int j, k = i0; k < i1; k++)
		{
			j = ja[k];
			rightVector[j] -= upperTringMat[k] * vectorX[i];
		}
	}
}

void SLAE::SolveBackwardLU(double* upperTringMat, double* rightVector, double* vectorX) {
	for (int i0, i1, i = n - 1; i >= 0; i--)
	{
		i0 = ia[i];
		i1 = ia[i + 1];
		vectorX[i] = rightVector[i];
		for (int j, k = i0; k < i1; k++)
		{
			j = ja[k];
			rightVector[j] -= upperTringMat[k] * vectorX[i];
		}
	}
}


void SLAE::MatrixUVectorMultiplicationLU(double* upperTringMat, double* vectorMult, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorMult[i];
		for (int j, k = ia[i]; k < ia[i + 1]; k++)
		{
			j = ja[k];
			vectorOut[j] += upperTringMat[k] * vectorMult[i];
		}
	}
}

void SLAE::MatrixUVectorMultiplicationLU(double* upperTringMat, double* diag, double* vectorMult, double* vectorOut)
{
	for (int i = 0; i < n; i++)
	{
		vectorOut[i] = vectorMult[i]*diag[i];
		for (int j, k = ia[i]; k < ia[i + 1]; k++)
		{
			j = ja[k];
			vectorOut[j] += upperTringMat[k] * vectorMult[i];
		}
	}
}

void SLAE::CalculateZ_LU(double* vectorOut)
{
	SolveBackwardLU(auLU, tmp1, tmp1);
	MatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(alLU, diLU, x0, tmp1);
	SolveBackwardLU(alLU, diLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, x0, vectorOut);
}

void SLAE::CalculateZ_LUaster(double* vectorOut)
{
	SolveBackwardLU(auLU, diLU, tmp1, tmp1);
	MatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(alLU, x0, tmp1);
	SolveBackwardLU(alLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, diLU, x0, vectorOut);
}

void SLAE::CalculateZ_LUsq(double* vectorOut)
{
	SolveBackwardLU(auLU, diLU, tmp1, tmp1);
	MatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(alLU, diLU, x0, tmp1);
	SolveBackwardLU(alLU, diLU, tmp1, tmp1);
	TransposedMatrixVectorMultiplication(tmp1, x0);
	SolveForwardLU(auLU, diLU, x0, vectorOut);
}

double SLAE::CalculateRelativeDiscrepancyWithR(double norm)
{
	return VectorNorm(r) / norm;
}

double SLAE::CalculateRelativeDiscrepancy(double norm)
{
	MatrixVectorMultiplication(x, tmp1);
	VectorSubtract(b, tmp1, tmp1);
	return VectorNorm(tmp1) / norm;
}

void SLAE::AllocateMemory()
{
	al = new double[nProfile]();
	au = new double[nProfile]();
	di = new double[n]();
	ja = new int[nProfile]();
	b = new double[n]();
	x = new double[n]();
	x0 = new double[n]();
	r = new double[n]();
	z = new double[n]();
	tmp1 = new double[n]();

	dP = new double[n]();

	xtrue = new double[n]();

	alLU = new double[nProfile]();
	auLU = new double[nProfile]();
	diLU = new double[n]();
}

void SLAE::ClearMemory()
{
	delete[] al;
	delete[] au;
	delete[] di;
	delete[] ja;
	delete[] b;
	delete[] x;
	delete[] x0;
	delete[] z;
	delete[] tmp1;

}