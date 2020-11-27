#include "System.h"
#include <iostream>
#include "Eigen/core"

using namespace std;

int main(void)
{
	System imusystem;
	
	cout << "AHRSEKF2 running " << endl;

	imusystem.RunEKF2();
}