#include <iostream>
#include "stdafx.h"
#include "FEM.h"

int main()
{
    FEM system;

    system.Input();
    system.Solve();
}