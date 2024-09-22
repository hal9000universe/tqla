//
// Created by Benjamin Schoofs on 11.09.24.
//

#ifndef MATRICOM_PERFORMANCE_H
#define MATRICOM_PERFORMANCE_H

#include <chrono>
#include "matrix.h"
#include "lossy_matrix.h"

using namespace std;
using namespace std::chrono;

void timeMatInversion();
void timeSparseAddition(int n);
void timeSparseMultiplication(int n);

#endif //MATRICOM_PERFORMANCE_H
