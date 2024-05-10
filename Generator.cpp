#include "stdafx.h"
#include "FEM.h"

void GenerateGrid(int factor)
{
    vector<Triangle> tris;
    vector<Vertex> verts;
	vector<vector<int>> bounds;
	double ta, tb, tStep;

	// Считываю узлы
	FILE* file;
	fopen_s(&file, "GenData/VerticesIN.txt", "r");

	fscanf_s(file, "%lf %lf %lf", &ta, &tStep, &tb);
	int num;
	fscanf_s(file, "%d", &num);
	Vertex tempVert{};

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%lf %lf", &tempVert.x, &tempVert.y);
		verts.push_back(tempVert);
	}
	fclose(file);
	// Считываю КЭ
	fopen_s(&file, "GenData/TrisIN.txt", "r");

	fscanf_s(file, "%d", &num);
	Triangle tempTri{};

	for (int i = 0; i < num; i++)
	{
		fscanf_s(file, "%d %d %d", &tempTri.vert1, &tempTri.vert2, &tempTri.vert3);
		tempTri.vert1--;
		tempTri.vert2--;
		tempTri.vert3--;
		tris.push_back(tempTri);
	}
	fclose(file);

	for (int f = 1; f <= factor; f++)
	{
		int end = tris.size();
		for (int i = 0; i < end; i++)
		{
			int newVInd[3]{};
			Triangle curTri = tris[0];
			tris.erase(tris.begin());

			Vertex newV[3]{};
			newV[0].x = (verts[curTri.vert1].x + verts[curTri.vert2].x) / 2.0;
			newV[0].y = (verts[curTri.vert1].y + verts[curTri.vert2].y) / 2.0;

			newV[1].x = (verts[curTri.vert2].x + verts[curTri.vert3].x) / 2.0;
			newV[1].y = (verts[curTri.vert2].y + verts[curTri.vert3].y) / 2.0;

			newV[2].x = (verts[curTri.vert3].x + verts[curTri.vert1].x) / 2.0;
			newV[2].y = (verts[curTri.vert3].y + verts[curTri.vert1].y) / 2.0;

			for (int v = 0; v < 3; v++)
			{
				bool found = false;
				for (int k = 0; k < verts.size(); k++)
				{
					if (abs(verts[k].x - newV[v].x) <= 1e-15 && abs(verts[k].y - newV[v].y) <= 1e-15) 
					{
						found = true;
						newVInd[v] = k;
					}
				}
				if (!found) {
					verts.push_back(newV[v]);
					newVInd[v] = verts.size() - 1;
				}
			}

			Triangle newT1, newT2, newT3, newT4;

			newT1.vert1 = newVInd[0];
			newT1.vert2 = newVInd[1];
			newT1.vert3 = newVInd[2];

			newT2.vert1 = curTri.vert1;
			newT2.vert2 = newVInd[0];
			newT2.vert3 = newVInd[2];

			newT3.vert1 = curTri.vert2;
			newT3.vert2 = newVInd[1];
			newT3.vert3 = newVInd[0];

			newT4.vert1 = curTri.vert3;
			newT4.vert2 = newVInd[2];
			newT4.vert3 = newVInd[1];

			tris.push_back(newT1);
			tris.push_back(newT2);
			tris.push_back(newT3);
			tris.push_back(newT4);
		}
	}
	for (int i = 0; i < verts.size(); i++)
	{
		if (abs(verts[i].x - verts[0].x) < 1e-15) 
		{
			vector<int> temp;
			temp.push_back(i + 1);
			temp.push_back(1);
			bounds.push_back(temp);
		}
		else if (abs(verts[i].y - verts[0].y) < 1e-15) {
			vector<int> temp;
			temp.push_back(i + 1);
			temp.push_back(2);
			bounds.push_back(temp);
		}
		else if (abs(verts[i].x - verts[1].x) < 1e-15) {
			vector<int> temp;
			temp.push_back(i + 1);
			temp.push_back(3);
			bounds.push_back(temp);
		}
		else if (abs(verts[i].y - verts[3].y) < 1e-15) {
			vector<int> temp;
			temp.push_back(i + 1);
			temp.push_back(4);
			bounds.push_back(temp);
		}
	}

	FILE* fileVerts;
	fopen_s(&fileVerts, "Vertices.txt", "w");
	fprintf_s(fileVerts, "%lf %lf %lf\n", ta, tStep / (pow(2, factor)), tb);
	fprintf_s(fileVerts, "%d\n", verts.size());
	for (int i = 0; i < verts.size(); i++)
		fprintf_s(fileVerts, "%.14lf %.14lf\n", verts[i].x, verts[i].y);
	fclose(fileVerts);

	FILE* fileTris;
	fopen_s(&fileTris, "Triangles.txt", "w");
	fprintf_s(fileTris, "%d\n", tris.size());
	for (int i = 0; i < tris.size(); i++)
		fprintf_s(fileTris, "%d %d %d 1\n", tris[i].vert1 + 1, tris[i].vert2 + 1, tris[i].vert3 + 1);
	fclose(fileTris);

	FILE* fileBound;
	fopen_s(&fileBound, "BoundaryConditions.txt", "w");
	fprintf_s(fileBound, "%d\n", bounds.size());
	for (int i = 0; i < bounds.size(); i++)
		fprintf_s(fileBound, "1 %d 1\n", bounds[i][0]);
	fclose(fileBound);
}

double AverageRate(vector<double> q0, vector<double> q1)
{
	double sum = 0;
	for (int i = 0, j = 6; i < q0.size(); i++, j += 4)
	{
		sum += log2(q0[i] / q1[j]);
	}
	return sum / q0.size();
}
