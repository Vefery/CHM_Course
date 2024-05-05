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
	vector<Vertex> vertices; // Узлы системы
	vector<Triangle> tris; // Конечные элементы
	vector<FirstBoundaryCondition> firstBoundary; // Первые краевые условия
	vector<SecondBoundaryCondition> secondBoundary; // Вторые краевые условия
	vector<ThirdBoundaryCondition> thirdBoundary; // Третьи краевые условия
	vector<double> timeStamps; // Узлы по времени
	vector<double> q; // Вектор решения
	vector<double> q_1, q_2; // Предыдущие векторы q(j-1) and q(j-2)

	void Input(); // Ввод данных из файла
	void Solve(); // Общая функция запуска решения
	void PrintSolution();
private:
	int regionsNum, globalN; // Количество областей и узлов
	vector<int> ig, jg; // Глобальная матрица
	vector<double> gglA, gguA, diA, gglM, gguM, diM, gglG, gguG, diG, b;
	double G[3][3]{}; // Пустая матрица G
	double M[3][3]{}; // Пустая матрица M
	const double pureM[3][3] = { {2, 1, 1}, {1, 2, 1}, {1, 1, 2} }; // Шаблон матрицы M для возвращения ее в исходное состояние на каждой итерации
	double localB[3]{}; // Локальный вектор b

	double Lamda(int vert, int region); // Вычисление значения лямбда
	double Gamma(int vert, int region); // Вычисление значения гамма
	double Function(int vert, int region); // Вычисление значения функции f
	double Beta(int vert, int eqNum); // Вычисление значения бета
	double Ubeta(int vert, int eqNum); // Вычисление значения U бета для 3 краевого условия
	double Theta(int vert, int eqNum); // Вычисление значения тета для 2 краевого условия
	double Ug(int vert, int eqNum); // Вычисление значения в узле для 1 краевого условия
	double GetAverageLamda(Triangle tri); // Вычисление среднего лямбда на элементе
	double GetAverageGamma(Triangle tri); // Вычисление среднего гамма на элементе
	double DetD(Triangle tri); // Вычисление определителя D (удвоенной площади) элемента
	double Alpha(Triangle tri, int k, int i);  // Вычисление значения альфа для построения матрицы G
	double EdgeLength(int vert1, int vert2); // Вычисление длины ребра
	int IndexOfUnknown(Triangle tri, int i); // Получение глобального номера узла из локального у элемента
	void FormM(Triangle tri); // Формирование матрицы G
	void FormG(Triangle tri); // Формирование матрицы M
	void FormPortrait(); // Формирование портрета глобальной матрицы
	void ResolveBoundaries(); // Учет всех краевых условий
	void AddToGlobal(int i, int j, double add); // Добавление значения в глобальную матрицу
	void AllocateGlobalMatrices(); // Выделение памяти для глобальной матрицы
	void FormB(Triangle tri); // Формирование локального вектора b
};