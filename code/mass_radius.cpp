#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

using namespace std;

int main() {


double M;
double R;
double maxRadius = 4;
double step = 0.1;
double Rearth = 1;
double test = 0.93;

ofstream valuesFile;
valuesFile.open ("Mass_Radius.csv");

for(R = 0; R < maxRadius; R += 0.1) {
	if (R < 1.5) {
		M = 0.440*pow(R/Rearth,3) + 0.614*pow(R/Rearth,4);
		valuesFile << M << ";" << R << endl;
		cout << M << ' ' << R << endl;
	}else{
		M = 2.69*pow(R,0.93);
		valuesFile << M << ";" << R << endl;
		cout << M << ' ' << R << endl;
	}
}
valuesFile.close();

return 0;

}
