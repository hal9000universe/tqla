//
// Created by Benjamin Schoofs on 12.09.24.
//

#include "lossy_matrix.h"

int getRowIndex(tuple<int, int, double> entry) {
    return get<0>(entry);
}

int getColumnIndex(tuple<int, int, double> entry) {
    return get<1>(entry);
}

double getValue(tuple<int, int, double> entry) {
    return get<2>(entry);
}

int getEntryIndex(vector<sparseEntry> entries, int row, int column) {
    for (int i = 0; i < entries.size(); ++i) {
        if (entries[i].row == row and entries[i].column == column) {
            return i;
        }
    }
    return -1;
}

lossyMatrix::lossyMatrix(vector<sparseEntry> entries, tuple<int, int> shape) {
    mEntries = entries;
    mShape = shape;
}

tuple<int, int> lossyMatrix::getShape() {
    return mShape;
}

vector<sparseEntry> lossyMatrix::getEntries() {
    return mEntries;
}

int lossyMatrix::numRows() {
    return get<0>(mShape);
}

int lossyMatrix::numColumns() {
    return get<1>(mShape);
}

double lossyMatrix::getEntry(int rowIndex, int columnIndex) {
    for (sparseEntry entry : mEntries) {
        if (entry.row == rowIndex and entry.column == columnIndex) {
            return entry.value;
        }
    }
    return 0;
}

lossyMatrix lossyMatrix::zeroMatrix(int numRows, int numColumns) {
    return lossyMatrix({}, make_tuple(numRows, numColumns));
}

matrix lossyMatrix::toMatrix() {
    vector<double> values;
    values.reserve(this->numRows() * this->numColumns());
    for (int i = 0; i < this->numRows(); ++i) {
        for (int j = 0; j < this->numColumns(); ++j) {
            values.push_back(this->getEntry(i, j));
        }
    }
    return matrix(values, mShape);
}

lossyMatrix lossyFromMatrix(matrix mat) {
    vector<sparseEntry> entries;
    for (int i = 0; i < mat.numRows(); ++i) {
        for (int j = 0; j < mat.numColumns(); ++j) {
            if (!doubleEqual(0, mat.getEntry(i, j))) {
                entries.push_back(sparseEntry(i, j, mat.getEntry(i, j)));
            }
        }
    }
    return lossyMatrix(entries, mat.getShape());
}

sparseEntry addEntries(sparseEntry aEntry, sparseEntry bEntry) {
    assert(aEntry.row == bEntry.row and aEntry.column == bEntry.column);
    return sparseEntry(aEntry.row, aEntry.column, aEntry.value + bEntry.value);
}

bool operator==(lossyMatrix &a, lossyMatrix &b) {
    assert(a.numRows() == b.numRows() and a.numColumns() == b.numColumns());
    if (a.getEntries().size() == b.getEntries().size()) {
        for (int i = 0; i < a.getEntries().size(); ++i) {
            sparseEntry aEntry = a.getEntries()[i], bEntry = b.getEntries()[i];
            if (aEntry.row != bEntry.row or aEntry.column != bEntry.column
                or !doubleEqual(aEntry.value, bEntry.value)) {
                return false;
            }
        }
    } else {
        return false;
    }
    return true;
}

lossyMatrix operator+(lossyMatrix &a, lossyMatrix &b) {
    assert(a.numRows() == b.numRows() and a.numColumns() == b.numColumns());

    vector<sparseEntry> aEntries = a.getEntries();
    vector<sparseEntry> bEntries = b.getEntries();
    vector<sparseEntry> outEntries;
    outEntries.reserve(aEntries.size() + bEntries.size());

    for (sparseEntry aEntry : aEntries) {
        bool added = false;
        for (int i = 0; i < bEntries.size(); i++) {
            if (bEntries[i].row == aEntry.row and bEntries[i].column == aEntry.column) {
                outEntries.push_back(aEntry + bEntries[i]);
                bEntries.erase(bEntries.begin() + i);
                added = true;
                break;
            } else if (aEntry.row < bEntries[i].row) {
                break;
            }
        }
        if (!added) {
            outEntries.push_back(aEntry);
        }
    }

    for (sparseEntry bEntry : bEntries) {
        outEntries.push_back(bEntry);
    }
    outEntries.shrink_to_fit();
    return lossyMatrix(outEntries, a.getShape());
}

lossyMatrix operator*(lossyMatrix &a, lossyMatrix &b) {
    assert(a.numColumns() == b.numRows());

    tuple<int, int> outShape = {a.numRows(), b.numColumns()};
    vector<tuple<int, int, double>> outEntries;
    vector<lossyMatrix> mats;
    mats.reserve(a.getEntries().size() * b.numColumns());

    for (sparseEntry aEntry : a.getEntries()) {
        for (sparseEntry bEntry : b.getEntries()) {
            if (aEntry.column == bEntry.row) {
                mats.push_back(lossyMatrix({aEntry * bEntry},outShape));
            } else if (aEntry.column < bEntry.row) {
                break;
            }
        }
    }

    lossyMatrix sumMatrix = lossyMatrix::zeroMatrix(a.numRows(), b.numColumns());
    for (lossyMatrix mat : mats) {
        sumMatrix = mat + sumMatrix;
    }
    return sumMatrix;
}