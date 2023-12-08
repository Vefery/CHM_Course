#include "stdafx.h"
#include "FEM.h"
#include "slae.h"

double FEM::Lamda(int vert, int region)
{
	switch (region)
	{
	case 1:
		return 1;
	case 2:
		return 1;
	default:
		return NAN;
	}
}

double FEM::Gamma(int vert, int region)
{
	switch (region)
	{
	case 1:
		return 1.0/3.0;
	case 2:
		return 1;
	default:
		return NAN;
	}
}

double FEM::Function(int vert, int region)
{
	switch (region)
	{
	case 1:
		return 4.0 - vertices[vert].x;
	case 2:
		return 6.0 * vertices[vert].x - 15.0;
	default:
		return NAN;
	}
}

double FEM::Beta(int vert, int eqNum)
{
	switch (eqNum)
	{
	case 1:
		return 1.0/5.0;
	case 2:
		return 1.0/6.0;
	default:
		return NAN;
	}
}

double FEM::Ubeta(int vert, int eqNum)
{
	switch (eqNum)
	{
	case 1:
		return 6.0 * vertices[vert].x;
	case 2:
		return -3.0 * vertices[vert].x;
	default:
		return NAN;
	}
}

double FEM::Theta(int vert, int eqNum)
{
	switch (eqNum)
	{
	case 1:
		return 0;
	case 2:
		return -6;
	case 3:
		return 0;
	default:
		return NAN;
	}
}

double FEM::Ug(int vert, int eqNum)
{
	switch (eqNum)
	{
	case 1:
		return 2.0 * vertices[vert].x + 5.0;
	default:
		return NAN;
	}
}

void FEM::Input()
{
	// —читывание узлов
	FILE* file;
	fopen_s(&file, "Vertices.txt", "r");

	int num;
	fscanf_s(file, "%d", &num);
	Vertex tempVert{};

	globalN = num;

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%lf %lf", &tempVert.x, &tempVert.y);
		vertices.push_back(tempVert);
	}
	fclose(file);

	// —читывание треугольников
	fopen_s(&file, "Triangles.txt", "r");

	fscanf_s(file, "%d", &num);
	Triangle tempTri{};

	regionsNum = num;

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%d %d %d %d", &tempTri.vert1, &tempTri.vert2, &tempTri.vert3, &tempTri.region);
		tempTri.vert1--;
		tempTri.vert2--;
		tempTri.vert3--;
		tris.push_back(tempTri);
	}
	fclose(file);

	// —читывание краевых
	fopen_s(&file, "BoundaryConditions.txt", "r");

	fscanf_s(file, "%d", &num);
	FirstBoundaryCondition firstBoundTemp{};
	SecondBoundaryCondition secondBoundTemp{};
	ThirdBoundaryCondition thirdBoundTemp{};

	for (int i = 0; i < num; i++)
	{
		int type;
		fscanf_s(file, "%d", &type);
		switch (type)
		{
		case 1:
			fscanf_s(file, "%d %d", &firstBoundTemp.vert, &firstBoundTemp.equationNum);
			firstBoundTemp.vert--;
			firstBoundary.push_back(firstBoundTemp);
			break;
		case 2:
			fscanf_s(file, "%d %d %d", &secondBoundTemp.vert1, &secondBoundTemp.vert2, &secondBoundTemp.equationNum);
			secondBoundTemp.vert1--;
			secondBoundTemp.vert2--;
			secondBoundary.push_back(secondBoundTemp);
			break;
		case 3:
			fscanf_s(file, "%d %d %d %d", &thirdBoundTemp.vert1, &thirdBoundTemp.vert2, &thirdBoundTemp.betaEquationNum, &thirdBoundTemp.UbetaEquationNum);
			thirdBoundTemp.vert1--;
			thirdBoundTemp.vert2--;
			thirdBoundary.push_back(thirdBoundTemp);
			break;
		default:
			break;
		}
	}
	fclose(file);
}

void FEM::Solve()
{
	FormPortrait();
	AllocateGlobalMatrix();
	for (int i = 0; i < regionsNum; i++)
	{
		Triangle currTri = tris[i];
		FormG(currTri);
		FormM(currTri);
		FormB(currTri);
		int globalBasis[3] = { currTri.vert1, currTri.vert2, currTri.vert3 };
		for (int j = 0; j < 3; j++)
		{
			b[globalBasis[j]] += localB[j];
			for (int k = 0; k < 3; k++)
				AddToGlobal(globalBasis[j], globalBasis[k], G[j][k] + M[j][k]);
		}
	}
	ResolveBoundaries();

	SLAE slae;
	slae.Input(globalN, 100000, 1e-15, ig, jg, ggl, ggu, di, b);
	slae.OutputDense();
	slae.MethodOfConjugateGradientsForNonSymMatrixWithDiagP();
	q = slae.x;
}

double FEM::GetAverageLamda(Triangle tri)
{
	return (Lamda(tri.vert1, tri.region) + Lamda(tri.vert2, tri.region) + Lamda(tri.vert3, tri.region)) / 3.0;
}

double FEM::GetAverageGamma(Triangle tri)
{
	return (Gamma(tri.vert1, tri.region) + Gamma(tri.vert2, tri.region) + Gamma(tri.vert3, tri.region)) / 3.0;
}

double FEM::DetD(Triangle tri)
{
	return (vertices[tri.vert2].x - vertices[tri.vert1].x) * (vertices[tri.vert3].y - vertices[tri.vert1].y) - (vertices[tri.vert3].x - vertices[tri.vert1].x) * (vertices[tri.vert2].y - vertices[tri.vert1].y);
}

double FEM::Alpha(Triangle tri, int k, int i)
{

	if (k == 1) {
		switch (i)
		{
		case 0:
			return vertices[tri.vert2].y - vertices[tri.vert3].y;
		case 1:
			return vertices[tri.vert3].y - vertices[tri.vert1].y;
		case 2:
			return vertices[tri.vert1].y - vertices[tri.vert2].y;
		default:
			return NAN;
		}
	}
	else if (k == 2) {
		switch (i)
		{
		case 0:
			return vertices[tri.vert3].x - vertices[tri.vert2].x;
		case 1:
			return vertices[tri.vert1].x - vertices[tri.vert3].x;
		case 2:
			return vertices[tri.vert2].x - vertices[tri.vert1].x;
		default:
			return NAN;
		}
	}
	else
		return NAN;
}

double FEM::EdgeLength(int vert1, int vert2)
{
	double x = vertices[vert1].x - vertices[vert2].x;
	double y = vertices[vert1].y - vertices[vert2].y;
	return sqrt(x*x + y*y);
}

int FEM::IndexOfUnknown(Triangle tri, int i)
{
	switch (i)
	{
	case 0:
		return tri.vert1;
	case 1:
		return tri.vert2;
	case 2:
		return tri.vert3;
	default:
		return NAN;
	}
}

void FEM::FormM(Triangle tri)
{
	M[0][0] = M[1][1] = M[2][2] = 2;
	M[0][1] = M[0][2] = M[1][0] = M[2][0] = M[1][2] = M[2][1] = 1;
	double temp = (GetAverageGamma(tri) * fabs(DetD(tri))) / 24.0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			M[i][j] *= temp;
	}
}

void FEM::FormG(Triangle tri)
{
	double averageLamda = GetAverageLamda(tri);
	double temp = averageLamda / (fabs(DetD(tri)) * 2.0);

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			G[i][j] = temp * (Alpha(tri, 1, i) * Alpha(tri, 1, j) + Alpha(tri, 2, i) * Alpha(tri, 2, j));
	}
}

void FEM::FormPortrait()
{
	ig = new int[globalN + 1];
	int* list[2]{};
	list[0] = new int[2 * globalN * (globalN - 2)];
	list[1] = new int[2 * globalN * (globalN - 2)];
	int* listbeg = new int[globalN];

	int listsize = -1;
	for (int i = 0; i < globalN; i++)
		listbeg[i] = -1;
	for (int ielem = 0; ielem < regionsNum; ielem++)
	{
		for (int i = 0; i < 3; i++)
		{
			int k = IndexOfUnknown(tris[ielem], i);
			for (int j = i + 1; j < 3; j++)
			{
				int ind1 = k;
				int ind2 = IndexOfUnknown(tris[ielem], j);
				if (ind2 < ind1) {
					ind1 = ind2;
					ind2 = k;
				}
				int iaddr = listbeg[ind2];
				if (iaddr == -1) {
					listsize++;
					listbeg[ind2] = listsize;
					list[0][listsize] = ind1;
					list[1][listsize] = -1;
				}
				else {
					while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
						iaddr = list[1][iaddr];

					if (list[0][iaddr] > ind1) {
						listsize++;
						list[0][listsize] = list[0][iaddr];
						list[1][listsize] = list[1][iaddr];
						list[0][iaddr] = ind1;
						list[1][iaddr] = listsize;
					}
					else {
						if (list[0][iaddr] < ind1) {
							listsize++;
							list[1][iaddr] = listsize;
							list[0][listsize] = ind1;
							list[1][listsize] = -1;
						}
					}
				}
			}
		}
	}

	jg = new int[listsize + 1];
	ig[0] = 0;
	for (int i = 0; i < globalN; i++)
	{
		ig[i + 1] = ig[i];
		int iaddr = listbeg[i];
		while (iaddr != -1) {
			jg[ig[i + 1]] = list[0][iaddr];
			ig[i + 1]++;
			iaddr = list[1][iaddr];
		}
	}
}

void FEM::ResolveBoundaries()
{
	// ”чет 3 краевых
	for (int i = 0; i < thirdBoundary.size(); i++)
	{
		ThirdBoundaryCondition temp = thirdBoundary[i];
		double localA[2][2] = { {2.0, 1.0}, {1.0, 2.0} };
		double factor = (((Beta(temp.vert1, temp.betaEquationNum) + Beta(temp.vert2, temp.betaEquationNum)) / 2.0) * EdgeLength(temp.vert1, temp.vert2)) / 6.0;
		b[temp.vert1] += factor * (2 * Ubeta(temp.vert1, temp.UbetaEquationNum) + Ubeta(temp.vert2, temp.UbetaEquationNum));
		b[temp.vert2] += factor * (Ubeta(temp.vert1, temp.UbetaEquationNum) + 2 * Ubeta(temp.vert2, temp.UbetaEquationNum));

		int globalBasis[2] = { temp.vert1, temp.vert2 };
		for (int i = 0; i < 2; i++)
		{
			for (int j = 0; j < 2; j++)
				AddToGlobal(globalBasis[i], globalBasis[j], localA[i][j] * factor);
		}
	}
	// ”чет 2 краевых
	for (int i = 0; i < secondBoundary.size(); i++)
	{
		SecondBoundaryCondition temp = secondBoundary[i];
		double factor = EdgeLength(temp.vert1, temp.vert2) / 6.0;
		double t1 = factor * (2 * Theta(temp.vert1, temp.equationNum) + Theta(temp.vert2, temp.equationNum));
		double t2 = factor * (Theta(temp.vert1, temp.equationNum) + 2 * Theta(temp.vert2, temp.equationNum));
		b[temp.vert1] += factor * (2 * Theta(temp.vert1, temp.equationNum) + Theta(temp.vert2, temp.equationNum));
		b[temp.vert2] += factor * (Theta(temp.vert1, temp.equationNum) + 2 * Theta(temp.vert2, temp.equationNum));
	}
	// ”чет 1 краевых
	for (int i = 0; i < firstBoundary.size(); i++)
	{
		FirstBoundaryCondition temp = firstBoundary[i];
		b[temp.vert] = Ug(temp.vert, temp.equationNum);

		di[temp.vert] = 1;
		for (int k = ig[temp.vert]; k < ig[temp.vert + 1]; k++) {
			ggl[k] = 0;
		}
		for (int k = 0; k < globalN; k++)
		{
			int ind;
			bool exist = false;
			for (ind = ig[k]; ind < ig[k + 1]; ind++)
			{
				if (jg[ind] == temp.vert) {
					exist = true;
					break;
				}

			}
			if (exist)
				ggu[ind] = 0;
		}
	}
}

void FEM::AddToGlobal(int i, int j, double add)
{
	if (i == j)
		di[i] += add;
	else if (i < j) {
		int ind;
		for (ind = ig[j]; ind < ig[j + 1]; ind++)
		{
			if (jg[ind] == i)
				break;
		}
		ggu[ind] += add;
	}
	else {
		int ind;
		for (ind = ig[i]; ind < ig[i + 1]; ind++)
		{
			if (jg[ind] == j)
				break;
		}
		ggl[ind] += add;
	}
}

void FEM::AllocateGlobalMatrix()
{
	di = new double[globalN]();
	b = new double[globalN]();

	ggl = new double[ig[globalN] - ig[0]]();
	ggu = new double[ig[globalN] - ig[0]]();
}

void FEM::FormB(Triangle tri)
{
	double factor = fabs(DetD(tri)) / 24.0;

	localB[0] = factor * (2.0 * Function(tri.vert1, tri.region) + Function(tri.vert2, tri.region) + Function(tri.vert3, tri.region));
	localB[1] = factor * (Function(tri.vert1, tri.region) + 2.0 * Function(tri.vert2, tri.region) + Function(tri.vert3, tri.region));
	localB[2] = factor * (Function(tri.vert1, tri.region) + Function(tri.vert2, tri.region) + 2.0 * Function(tri.vert3, tri.region));
}