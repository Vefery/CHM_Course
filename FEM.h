#pragma once
#include "stdafx.h"

using namespace std;

struct Vertex
{
	double x;
	double y;
};

struct Triangle
{
	int vert1, vert2, vert3;
	int region;
};

struct FirstBoundaryCondition 
{
	int vert;
	int equationNum;
};

struct SecondBoundaryCondition
{
	int vert1, vert2;
	int equationNum;
};

struct ThirdBoundaryCondition
{
	int vert1, vert2;
	int betaEquationNum, UbetaEquationNum;
};

class FEM
{
public:
	vector<Vertex> vertices;
	vector<Triangle> tris;
	vector<FirstBoundaryCondition> firstBoundary;
	vector<SecondBoundaryCondition> secondBoundary;
	vector<ThirdBoundaryCondition> thirdBoundary;
	double* q;

	void Input();
	void Solve();
private:
	int regionsNum, globalN;
	int* ig, * jg;
	double* ggl, * ggu, * di, * b;
	double G[3][3]{};
	double M[3][3]{};
	const double pureM[3][3] = { {2, 1, 1}, {1, 2, 1}, {1, 1, 2} };
	double localB[3]{};

	double Lamda(int vert, int region);
	double Gamma(int vert, int region);
	double Function(int vert, int region);
	double Beta(int vert, int eqNum);
	double Ubeta(int vert, int eqNum);
	double Theta(int vert, int eqNum);
	double Ug(int vert, int eqNum);
	double GetAverageLamda(Triangle tri);
	double GetAverageGamma(Triangle tri);
	double DetD(Triangle tri);
	double Alpha(Triangle tri, int k, int i);
	double EdgeLength(int vert1, int vert2);
	int IndexOfUnknown(Triangle tri, int i);
	void FormM(Triangle tri);
	void FormG(Triangle tri);
	void FormPortrait();
	void ResolveBoundaries();
	void AddToGlobal(int i, int j, double add);
	void AllocateGlobalMatrix();
	void FormB(Triangle tri);
};