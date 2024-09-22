//
// Created by Benjamin Schoofs on 11.09.24.
//

#include "matrix.h"

matrix::matrix(vector<double> matrixValues, tuple<int, int> matrixShape) {
    assert(matrixValues.size() == get<0>(matrixShape) * get<1>(matrixShape));
    mValues = std::move(matrixValues);
    mShape = matrixShape;
}

tuple<int, int> matrix::getShape() {
    return mShape;
}

int matrix::numRows() {
    return get<0>(mShape);
}

int matrix::numColumns() {
    return get<1>(mShape);
}

double matrix::getEntry(int rowIndex, int columnIndex) {
    int index = this->numRows() * rowIndex + columnIndex;
    return mValues[index];
}

vector<double> matrix::getEntries() {
    return mValues;
}

int matrix::length() {
    return this->numRows() * this->numColumns();
}

matrix matrix::zero(int dim) {
    vector<double> zeroVector;
    zeroVector.reserve(dim * dim);
    tuple<double, double> zeroShape = {dim, dim};
    for (int i = 0; i < dim * dim; ++i) {
        zeroVector.push_back(0.);
    }
    return matrix{zeroVector, zeroShape};
}

matrix matrix::identity(int dim) {
    vector<double> identityVector;
    identityVector.reserve(dim * dim);
    tuple<double, double> identityShape = {dim, dim};
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            if (i == j) {
                identityVector.push_back(1.);
            } else {
                identityVector.push_back(0.);
            }
        }
    }
    return matrix{identityVector, identityShape};
}

matrix matrix::copy() {
    vector<double> copyVector = mValues;
    return matrix{copyVector, mShape};
}

bool matrix::quadratic() {
    return (this->numRows() == this->numColumns());
}

bool matrix::invertible() {
    if (abs(this->determinant()) < 1e-10) {
        return false;
    } else {
        return true;
    }
};

bool matrix::upperTriangular() {
    for (int i = 0; i < this->numRows(); ++i) {
        for (int j = 0; j < this->numColumns(); ++j) {
            if (i > j and abs(this->getEntry(i, j)) > 1e-10) {
                return false;
            }
        }
    }
    return true;
}

bool matrix::lowerTriangular() {
    for (int i = 0; i < this->numRows(); ++i) {
        for (int j = 0; j < this->numColumns(); ++j) {
            if (i < j and abs(this->getEntry(i, j)) > 1e-10) {
                return false;
            }
        }
    }
    return true;
}

void matrix::swapRows(int a, int b) {
    for (int j = 0; j < this->numColumns(); ++j) {
        double aEntry = this->getEntry(a, j);
        double bEntry = this->getEntry(b, j);
        mValues[this->numRows() * a + j] = bEntry;
        mValues[this->numRows() * b + j] = aEntry;
    }
}

void matrix::rowAdd(int addedIndex, int changedIndex, double factor) {
    for (int j = 0; j < this->numColumns(); ++j) {
        mValues[this->numRows() * changedIndex + j] =
                this->getEntry(changedIndex, j)
                + factor * this->getEntry(addedIndex, j);
    }
}

void matrix::columnAdd(int addedIndex, int changedIndex, double factor) {
    for (int i = 0; i < this->numRows(); ++i) {
        mValues[this->numRows() * i + changedIndex] =
                this->getEntry(i, changedIndex)
                + factor * this->getEntry(i, addedIndex);
    }
}

void matrix::rowMul(int changedIndex, double factor) {
    for (int j = 0; j < this->numColumns(); ++j) {
        mValues[this->numRows() * changedIndex + j] = factor * this->getEntry(changedIndex, j);
    }
}

void matrix::columnMul(int changedIndex, double factor) {
    for (int i = 0; i < this->numRows(); ++i) {
        mValues[this->numRows() * i + changedIndex] = factor * this->getEntry(i, changedIndex);
    }
}

matrix matrix::transpose() {
    vector<double> transposeVector;
    transposeVector.reserve(this->length());
    tuple<double, double> transposeShape = {this->numColumns(), this->numRows()};
    for (int i = 0; i < this->numRows(); ++i) {
        for (int j = 0; j < this->numColumns(); ++j) {
            transposeVector.push_back(this->getEntry(j, i));
        }
    }
    matrix transpose(transposeVector, transposeShape);
    return transpose;
}

matrix matrix::minor(int rowIndex, int columnIndex) {
    int minorLength = (this->numRows() - 1) * (this->numColumns() - 1);
    vector<double> minorVector;
    minorVector.reserve(minorLength);
    tuple<int, int> minorShape = {this->numRows() - 1, this->numColumns() - 1};
    for (int i = 0; i < this->numRows(); ++i) {
        for (int j = 0; j < this->numColumns(); ++j) {
            if (i == rowIndex or j == columnIndex) {
                continue;
            } else {
                minorVector.push_back(this->getEntry(i, j));
            }
        }
    }
    return matrix{minorVector, minorShape};
}

double matrix::determinant() {
    assert(this->quadratic());

    if (this->upperTriangular() or this->lowerTriangular()) {
        double det = 1.;
        for (int i = 0; i < this->numRows(); ++i) {
            det *= this->getEntry(i, i);
        }
        return det;
    } else if (this->numRows() == 2) {
        return this->getEntry(0, 0) * this->getEntry(1, 1)
               - this->getEntry(1, 0) * this->getEntry(0, 1);
    } else {
        double det = 0.;
        for (int i = 0; i < this->numRows(); ++i) {
            matrix minor = this->minor(i, 0);
            det += pow(-1., i) * this->getEntry(i, 0) * minor.determinant();
        }
        return det;
    }
}

matrix matrix::gaussianInverse() {
    matrix identity = matrix::identity(this->numRows());
    matrix inverse = matrix::identity(this->numRows());
    matrix helper = this->copy();

    // forward
    if (!this->upperTriangular()) {
        for (int j = 0; j < this->numColumns(); ++j) {
            // create non-zero element
            if (abs(helper.getEntry(j, j)) < 1e-10) {
                for (int i = j + 1; i < this->numRows(); ++i) {
                    if (abs(this->getEntry(i, j)) > 1e-10) {
                        helper.swapRows(j, i);
                        inverse.swapRows(j, i);
                        break;
                    }
                }

                if (abs(helper.getEntry(j, j)) < 1e-10) {
                    continue;
                }
            }

            // normalize j-th row
            if (abs(helper.getEntry(j, j)) > 1e-10) {
                double normalizer = 1. / helper.getEntry(j, j);
                helper.rowMul(j, normalizer);
                inverse.rowMul(j, normalizer);
            } else {
                return matrix::zero(1);
            }

            // eliminate lower rows
            for (int i = j + 1; i < this->numRows(); ++i) {
                double factor = -1. / helper.getEntry(i, j);
                helper.rowAdd(j, i, factor);
                inverse.rowAdd(j, i, factor);
            }
        }
    }

    // backward
    if (!this->upperTriangular()) {
        for (int j = this->numColumns(); j >= 0; --j) {
            if (abs(helper.getEntry(j, j)) > 1e-10) {
                double normalizer = 1. / helper.getEntry(j, j);
                helper.rowMul(j, normalizer);
                inverse.rowMul(j, normalizer);
            } else {
                return matrix::zero(0);
            }

            // eliminate upper rows
            for (int i = 0; i <= j; ++i) {
                double factor = -1. / helper.getEntry(i, j);
                helper.rowAdd(j, i, factor);
                inverse.rowAdd(j, i, factor);
            }
        }
    }

    return inverse;
}

matrix matrix::cramerInverse() {
    double det = this->determinant();
    if (abs(det) > 1e-10) {
        vector<double> inverseVector;
        inverseVector.reserve(this->length());
        tuple<double, double> inverse_shape = mShape;
        for (int i = 0; i < this->numRows(); ++i) {
            for (int j = 0; j < this->numColumns(); ++j) {
                inverseVector.push_back(pow(-1, i + j) *
                                        this->minor(i, j).determinant() / det);

            }
        }
        matrix preCramerInverse(inverseVector, inverse_shape);
        return preCramerInverse.transpose();
    } else {
        return matrix::zero(0);
    }
}

void matrix::print() {
    cout << "[";
    for (int i = 0; i < this->numRows(); ++i) {
        if (i == 0) {
            cout << "[";
        } else {
            cout << " [";
        }
        for (int j = 0; j < this->numColumns(); ++j) {
            cout << this->getEntry(i, j);
            if (j == this->numColumns() - 1) {
                if (i == this->numRows() - 1) {
                    cout << "]]" << endl;
                } else {
                    cout << "]" << endl;
                }
            } else {
                cout << ", ";
            }
        }
    }
}

bool operator==(matrix &a, matrix &b) {
    assert(a.numRows() == b.numRows() and a.numColumns() == b.numColumns());
    for (int i = 0; i < a.numRows(); ++i) {
        for (int j = 0; j < a.numColumns(); ++j) {
            if (!doubleEqual(a.getEntry(i, j), b.getEntry(i, j))) {
                return false;
            }
        }
    }
    return true;
}

matrix operator+(matrix &a, matrix &b) {
    assert(a.numRows() == b.numRows() and a.numColumns() == b.numColumns());

    tuple<int, int> cShape = {a.numRows(), a.numColumns()};
    vector<double> cVector;
    cVector.reserve(a.length());

    for (int i = 0; i < a.numRows(); ++i) {
        for (int j = 0; j < a.numColumns(); ++j) {
            cVector.push_back(a.getEntry(i, j) + b.getEntry(i, j));
        }
    }
    return matrix{cVector, cShape};
}

matrix operator*(matrix &a, matrix &b) {
    assert(a.numColumns() == b.numRows());

    tuple<int, int> cShape = {a.numRows(), b.numColumns()};
    int cLength = a.numRows() * b.numColumns();
    vector<double> cVector;
    cVector.reserve(cLength);

    for (int i = 0; i < a.numRows(); ++i) {
        for (int j = 0; j < b.numColumns(); ++j) {
            double entry = 0.;
            for (int k = 0; k < a.numColumns(); ++k) {
                entry += a.getEntry(i, k) * b.getEntry(k, j);
            }
            cVector.push_back(entry);
        }
    }
    return matrix{cVector, cShape};
}