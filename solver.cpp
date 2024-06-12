//-static-libgcc -lgdi32 -lcomdlg32 -I"C:\Program Files (x86)\eigen-3.4.0"

//Frequently used
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

//Task specific
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
using namespace Eigen;

// Function to read the boundary values from the file
void readBoundaryValues(const string &filename, int &steps, vector<double> &boundaryValues) {
    ifstream inputFile(filename);
    if (!inputFile) {
        cerr << "Error opening file" << endl;
        exit(1);
    }

    inputFile >> steps;
    int totalBoundaryPoints = 4 * steps;
    boundaryValues.resize(totalBoundaryPoints);

    for (int i = 0; i < totalBoundaryPoints; ++i) {
        inputFile >> boundaryValues[i];
    }

    inputFile.close();
}//


int main(){
	string filename = "boundary_values.txt";
	int steps;
	vector<double> boundaryValues;
    readBoundaryValues(filename, steps, boundaryValues);
    
    //Size of the grid
    int n = steps + 1;
    int gridPoints = n * n;
    
    // Create the sparse matrix and the right-hand side vector
    SparseMatrix<double> A(gridPoints, gridPoints);
    VectorXd b = VectorXd::Zero(gridPoints);
    
	vector<Triplet<double>> coefficients;
	// Fill the sparse matrix and the RHS vector
	// Lambda function idx for indices
    auto idx = [&](int i, int j) { return i * n + j; };

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int index = idx(i, j);

            if (i == 0 || i == n - 1 || j == 0 || j == n - 1) {
                // Boundary points
                //b(index) = boundaryValues[i == 0 ? j : i == n - 1 ? n + j - 1 : j == 0 ? 3 * n - i - 1 : 2 * n + n - i - 1];
                b(index) = boundaryValues[
					j==0 ? i :
					i==n-1 ? (n-1) + j :
					j==n-1 ? 2*(n - 1) + (n-1 - i) :
					3*(n-1) + (n-1 - j)
					
					];
					
                	
                coefficients.push_back(Triplet<double>(index, index, 1.0));
            } else {
                // Interior points
                coefficients.push_back(Triplet<double>(index, idx(i - 1, j), 1.0));
                coefficients.push_back(Triplet<double>(index, idx(i + 1, j), 1.0));
                coefficients.push_back(Triplet<double>(index, idx(i, j - 1), 1.0));
                coefficients.push_back(Triplet<double>(index, idx(i, j + 1), 1.0));
                coefficients.push_back(Triplet<double>(index, index, -4.0));
            }
        }
    }

    // Assemble the sparse matrix
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    // Solve the linear system
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(A);
    if (solver.info() != Success) {
        cerr << "Decomposition failed" << endl;
        return -1;
    }

    VectorXd x = solver.solve(b);
    if (solver.info() != Success) {
        cerr << "Solving failed" << endl;
        return -1;
    }

    // Output the solution
//    for (int i = 0; i < n; ++i) {
//        for (int j = 0; j < n; ++j) {
//            cout << x[idx(i, j)] << " ";
//        }
//        cout << endl;
//    }

	for (int i = 0; i < n; ++i) {
    	for (int j = 0; j < n; ++j) {
        	double xi = static_cast<double>(i) / (n - 1);
            double yj = static_cast<double>(j) / (n - 1);
            std::cout << xi << " " << yj << " " << x[idx(i, j)] << std::endl;
    	}
	}

    return 0;
}
	
	
	