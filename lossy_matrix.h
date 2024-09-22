//
// Created by Benjamin Schoofs on 12.09.24.
//

#ifndef MATRICOM_LOSSY_MATRIX_H
#define MATRICOM_LOSSY_MATRIX_H

#include <iostream>
#include <vector>
#include <cassert>
#include <cmath>
#include <chrono>

#include "helper.h"
#include "matrix.h"
#include "sparse_entry.h"

using namespace std;
using namespace std::chrono;

int getRowIndex(tuple<int, int, double> &entry);
int getColumnIndex(tuple<int, int, double> &entry);
double getValue(tuple<int, int, double> entry);
int getEntryIndex(vector<sparseEntry> entries, int row, int column); // TODO: to optimize

class lossyMatrix {
    private:
        vector<sparseEntry> mEntries;
        tuple<int, int> mShape;

    public:
        lossyMatrix(vector<sparseEntry> entries, tuple<int, int> shape);
        tuple<int, int> getShape();
        vector<sparseEntry> getEntries();
        int numRows();
        int numColumns();
        double getEntry(int rowIndex, int columnIndex);
        static lossyMatrix zeroMatrix(int numRows, int numColumns);
        matrix toMatrix();
};

lossyMatrix lossyFromMatrix(matrix mat);
sparseEntry addEntries(sparseEntry aEntry, sparseEntry bEntry);

bool operator==(lossyMatrix &a, lossyMatrix &b);
lossyMatrix operator+(lossyMatrix &a, lossyMatrix &b);
lossyMatrix operator*(lossyMatrix &a, lossyMatrix &b);

#endif //MATRICOM_LOSSY_MATRIX_H
