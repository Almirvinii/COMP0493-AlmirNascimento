#include <bits/stdc++.h>
using namespace std;

// Method to calculate prefix sum till index idx
int sum(int idx, vector<int>& F)
{
    int runningSum = 0;
    // Summing up all the partial sums
    while (idx > 0) {
        runningSum += F[idx];
        int rightMostSetBit = (idx & (-idx));
        idx -= rightMostSetBit;
    }
    return runningSum;
}

// Method to update the array by adding X to index idx
void add(int idx, int X, vector<int>& F)
{
    while (idx < F.size()) {
        F[idx] += X;
        int rightMostSetBit = (idx & (-idx));
        idx += rightMostSetBit;
    }
}

// Method to calculate the prefix sum till ith index
int prefSum(int idx, vector<int> &F1, vector<int> &F2) {
    return sum(idx, F1) * idx - sum(idx, F2);
}

// Method to calculate sum of range [L, R]
int rangeQuery(int L, int R, vector<int>& F1, vector<int> &F2) {
    return prefSum(R, F1, F2) - prefSum(L - 1, F1, F2);
}

// Add X to all elements in range [l, r]
void rangeUpdate(int l, int r, int X, vector<int> &F1, vector<int> &F2) {
    add(l, X, F1);
    add(r + 1, -X, F1);
    add(l, X * (l - 1) , F2);
    add(r + 1, -(X * r), F2);
}