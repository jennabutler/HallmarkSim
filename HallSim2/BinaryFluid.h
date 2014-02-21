/*
 * BinaryFluid.h
 *
 *  Created on: Feb 13, 2012
 *      Author: jenna
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <sstream>

#ifndef BINARYFLUID_H_
#define BINARYFLUID_H_

class BinaryFluid {
public:
	BinaryFluid();
	BinaryFluid(int cellGridSize);
	void updateBinaryFluid();
	void updateBinaryFluidSmall();
	void parametercalc();
	void streamout();
	void equilibriumdist();
	void initialize();
	double getOxygenValue(int i, int j);
	std::string formatOxygenValue();
	void consumeOxygen(int i, int j, double amount);
	std::string outputOxygenGrid(int iter);
	void smallConsumeOxygen();
	void printOxygen();
	void zeroFf();
	void rescale();
	//virtual ~BinaryFluid();

};

#endif /* BINARYFLUID_H_ */
