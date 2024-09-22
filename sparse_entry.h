//
// Created by Benjamin Schoofs on 13.09.24.
//

#ifndef MATRICOM_SPARSE_ENTRY_H
#define MATRICOM_SPARSE_ENTRY_H

#include <cassert>

class sparseEntry {
    public:
    int row;
    int column;
    double value;
    sparseEntry(int r, int c, double v);
};

bool operator<(sparseEntry a, sparseEntry b);
sparseEntry operator+(sparseEntry a, sparseEntry b);
sparseEntry operator*(sparseEntry a, sparseEntry b);

#endif //MATRICOM_SPARSE_ENTRY_H
