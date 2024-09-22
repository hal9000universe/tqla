//
// Created by Benjamin Schoofs on 11.09.24.
//

#include "performance.h"

void timeMatInversion() {
    cout << "Timing Matrix Inversion: " << endl;
    vector<double> vec = {1., 2., 0., 6.,
                          2., 5., 17, 0.,
                          0., 0., 10., 3.,
                          0., 0., 3., 1.};
    tuple<int, int> shape = {4, 4};
    matrix mat(vec, shape);

    // time cramer inversion
    auto cStart = high_resolution_clock::now();
    matrix cramerInverse = mat.cramerInverse();
    auto cStop = high_resolution_clock::now();
    auto cDuration = duration_cast<microseconds>(cStop - cStart);
    cout << "Cramer Duration: " << cDuration.count() << " ms" << endl;

    // time safe gaussian inversion
    auto gStart = high_resolution_clock::now();
    matrix gaussianInverse = mat.gaussianInverse();
    auto gStop = high_resolution_clock::now();
    auto gDuration = duration_cast<microseconds>(gStop - gStart);
    cout << "Gauss Duration: " << gDuration.count() << " ms" << endl;
    cout << "-------------- Finished! --------------" << endl;
}

void timeSparseAddition(int n) {
    cout << "Timing Matrix Addition: " << endl;
    // dense
    matrix mat = matrix::identity(n);
    auto matrixStart = high_resolution_clock::now();
    matrix denseResult = mat + mat;
    auto matrixStop = high_resolution_clock::now();
    auto matrixDuration = duration_cast<microseconds>(matrixStop - matrixStart);
    cout << "Dense Duration: " << matrixDuration.count() << " ms" << endl;
    // lossy
    lossyMatrix lossyMat = lossyFromMatrix(mat);
    auto lossyStart = high_resolution_clock::now();
    lossyMatrix lossyResult = lossyMat + lossyMat;
    auto lossyStop = high_resolution_clock::now();
    auto lossyDuration = duration_cast<microseconds>(lossyStop - lossyStart);
    cout << "Lossy Duration: " << lossyDuration.count() << " ms" << endl;
    cout << "-------------- Finished! --------------" << endl;
}

void timeSparseMultiplication(int n) {
    cout << "Timing Matrix Multiplication: " << endl;
    // dense
    matrix mat = matrix::identity(n);
    auto matrixStart = high_resolution_clock::now();
    matrix denseResult = mat * mat;
    auto matrixStop = high_resolution_clock::now();
    auto matrixDuration = duration_cast<microseconds>(matrixStop - matrixStart);
    cout << "Dense Duration: " << matrixDuration.count() << " ms" << endl;
    // lossy
    lossyMatrix lossyMat = lossyFromMatrix(mat);
    auto lossyStart = high_resolution_clock::now();
    lossyMatrix lossyResult = lossyMat * lossyMat;
    auto lossyStop = high_resolution_clock::now();
    auto lossyDuration = duration_cast<microseconds>(lossyStop - lossyStart);
    cout << "Lossy Duration: " << lossyDuration.count() << " ms" << endl;
    cout << "-------------- Finished! --------------" << endl;
}