//
// Created by Benjamin Schoofs on 13.09.24.
//

#include "sparse_entry.h"

sparseEntry::sparseEntry(int r, int c, double v) {
    row = r;
    column = c;
    value = v;
}

bool operator<(sparseEntry a, sparseEntry b) {
    return (a.row < b.row or (a.row == b.row and a.column < b.column));
}

sparseEntry operator+(sparseEntry a, sparseEntry b) {
    assert(a.row == b.row and a.column == b.column);
    return {a.row, a.column, a.value + b.value};
}

sparseEntry operator*(sparseEntry a, sparseEntry b) {
    assert(a.column == b.row);
    return {a.row, b.column, a.value * b.value};
}