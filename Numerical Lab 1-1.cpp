// Numerical Lab 1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//Name: Kizitor Chukwuma
//Date: 9/12/2019

#include "pch.h"
#include <iostream>
#include <math.h>
using namespace std;

int rows = 0;
int cols = 0;

double* x = new double[rows];

void BackSub(double** matrix) //Answer for number 2
{
	for (int i = rows - 1; i >= 0; i--)                
	{                        
		x[i] = matrix[i][rows];                
		for (int j = i + 1; j < rows; j++)
		{
			if (j != i)
			{                                  
				x[i] = x[i] - matrix[i][j] * x[j];
			}
		}
		x[i] = x[i] / matrix[i][i];            
	}
}

double** NaiveGauss(double** matrix) //Answer for 1
{
	for (int i = 0; i < rows - 1; i++)
	{
		for (int k = i + 1; k < rows; k++)
		{
			double t = matrix[k][i] / matrix[i][i];
			for (int j = 0; j <= rows; j++)
				matrix[k][j] = matrix[k][j] - t * matrix[i][j];    
		}
	}
	BackSub(matrix);
	return matrix;
}



int main() //Answer for 3
{
	rows = 4;
	cols = 5;
	cout << rows;
	cout << cols;
	double** matrix = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		matrix[i] = new double[cols]; 
	}

	matrix[0][0] = 1;
	matrix[0][1] = 3;
	matrix[0][2] = 2;
	matrix[0][3] = 1;
	matrix[0][4] = -2;

	matrix[1][0] = 4;
	matrix[1][1] = 2;
	matrix[1][2] = 1;
	matrix[1][3] = 2;
	matrix[1][4] = 2;

	matrix[2][0] = 2;
	matrix[2][1] = 1;
	matrix[2][2] = 2;
	matrix[2][3] = 3;
	matrix[2][4] = 1;

	matrix[3][0] = 1;
	matrix[3][1] = 2;
	matrix[3][2] = 4;
	matrix[3][3] = 1;
	matrix[3][4] = -1;

	for (int i = 0; i < rows; i++)
	{
		cout << endl;
		for (int j = 0; j < cols; j++)
		{
			cout << matrix[i][j] << " ";
		}
	}

	matrix = NaiveGauss(matrix);

	cout << endl;

	for (int i = 0; i < rows; i++)
	{
		cout << endl;
		for (int j = 0; j < cols; j++)
		{
			cout << matrix[i][j] << " ";
		}
	}

	cout << endl;

	for (int i = 0; i < rows; i++) //Answer for part a
	{
		cout << endl;
		cout << x[i];
	}

	cout << endl;

	for (int i = 1; i < 16; i++) //Answer for part b
	{
		cout << pow(2.0, -4.0*i) << " ";
	}

	cout << endl;
	cout << "When the varaiable gets small the error becomes worse due to round off error";
}