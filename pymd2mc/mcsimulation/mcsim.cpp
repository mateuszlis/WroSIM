/*
 * mcsim.cpp
 *
 *  Created on: 23-04-2011
 *      Author: lisu
 */

#include <iostream>

#include <cstdlib>
#include <ctime>
#include <cmath>
#include "TriangularLattice.h"
#include <fstream>

using namespace std;

const int LATT_SIZE = 10000;
const int EQUIB_STEPS = 2000;
const double R = 1.986;
const int T = 310;

double calcEnergyDiff(int *latt, int lattSize, int rowSize, int pos1, int pos2,
		double omegaAB);
void metropolis(TriangularLattice *latt, double omegaAB, int stepSize,
		int steps, ostream &stream, ostream &neighOstream);
void printLatt(TriangularLattice *lattice, int lattSize, int rowSize,
		ostream &output);

int main() {
	srand((unsigned) time(0));
	float omega = -1000;
	int aLipidsNum = 2000;

	ofstream trajFile, neighHistFile;
	trajFile.open("traj_omega-1000.xyz");
	neighHistFile.open("neigh_hist_omega-1000.dat");

	TriangularLattice *lattice = new TriangularLattice(LATT_SIZE, 100,
			aLipidsNum);
	metropolis(lattice, omega, 10000, 10000, trajFile, neighHistFile);
	trajFile.close();
}

double prob(double dG) {
	return exp(-dG / (R * T));
}

double calcEnergyDiff(TriangularLattice *latt, int pos1, int pos2,
		double omegaAB) {
	int diff1 = (6 - latt->simNeighbCount(pos1)) + (6 - latt->simNeighbCount(
			pos2));
	latt->exchangeSites(pos1, pos2);
	int diff2 = (6 - latt->simNeighbCount(pos1)) + (6 - latt->simNeighbCount(
			pos2));
	latt->exchangeSites(pos1, pos2);
	return (diff2 - diff1) * omegaAB;
}
void metropolisStep(TriangularLattice *latt, double omegaAB) {
	int pos1 = rand() % latt->getLatticeSize();
	int pos2 = rand() % latt->getLatticeSize();
	if ((*latt)[pos1] != (*latt)[pos2]) {
		double p = prob(calcEnergyDiff(latt, pos1, pos2, omegaAB));
		double acceptance = rand() / (float(RAND_MAX) + 1);
		if (p >= 1 or p > (acceptance)) {
			//cout << "DONE MOVE " << pos1 << " " << pos2 << endl;
			latt->exchangeSites(pos1, pos2);
		}
	}
}

void createNeighHist(long long *histArr, TriangularLattice *latt) {
	for (int i = 0; i < latt->getLatticeSize(); i++) {
		if ((*latt)[i])
			histArr[latt->simNeighbCount(i)] += 1;
	}
}
void metropolis(TriangularLattice *latt, double omegaAB, int stepSize,
		int steps, ostream &stream, ostream &neighOstream) {
	long long neighHist[7] = { 0, 0, 0, 0, 0, 0, 0 };
	for (int i = 0; i < steps; i++) {
		for (int j = 0; j < stepSize; j++) {
			metropolisStep(latt, omegaAB);

		}

		if (i % 10 == 0) {
			stream << (*latt);
		}
		if (i > EQUIB_STEPS) {
			createNeighHist(neighHist, latt);
		}

	}
	for (int i = 0; i < 7; i++) {
		neighOstream << i << " " << neighHist[i] / (steps - EQUIB_STEPS)
				<< endl;
	}

}
