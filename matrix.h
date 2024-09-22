//
// Created by Benjamin Schoofs on 11.09.24.
//

#ifndef MATRICOM_MATRIX_H
#define MATRICOM_MATRIX_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>

#include "helper.h"

using namespace std;

class matrix {
private:
    vector<double> mValues;
    tuple<int, int> mShape;

public:
    matrix(vector<double> matrixValues, tuple<int, int> matrixShape);
    tuple<int, int> getShape();
    int numRows();
    int numColumns();
    double getEntry(int rowIndex, int columnIndex);
    vector<double> getEntries();
    int length();
    static matrix zero(int dim);
    static matrix identity(int dim);
    matrix copy();
    bool quadratic();
    bool invertible();
    bool upperTriangular();
    bool lowerTriangular();
    void swapRows(int a, int b);
    void rowAdd(int addedIndex, int changedIndex, double factor);
    void columnAdd(int addedIndex, int changedIndex, double factor);
    void rowMul(int changedIndex, double factor);
    void columnMul(int changedIndex, double factor);
    matrix transpose();
    matrix minor(int rowIndex, int columnIndex);
    double determinant();
    matrix gaussianInverse();
    matrix cramerInverse();
    void print();
};

bool operator==(matrix &a, matrix &b);
matrix operator+(matrix &a, matrix &b);
matrix operator*(matrix &a, matrix &b);

#endif //MATRICOM_MATRIX_H