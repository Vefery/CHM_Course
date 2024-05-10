#include <iostream>
#include "stdafx.h"
#include "FEM.h"
//#include "Generator.h"

int main()
{
    FEM fem, fem2;
    int startDiv = 1;
    vector<double> q0, q1;
    // На сетке N
    //GenerateGrid(startDiv);
    fem.Input();
    fem.Solve();
    //q0 = fem.errVec;

    //printf_s("\n\n");
    //// На сетке 2N
    //GenerateGrid(startDiv + 1);
    //fem2.Input();
    //fem2.Solve();
    //q1 = fem2.errVec;

    //printf_s("%lf", log2(q0[q0.size() - 1] / q1[q1.size() - 1]));
}