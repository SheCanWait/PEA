#pragma once
#include "FileReader.h"

class PEA_0 
{
public:
	static void calculateSolution(FileReader* fileReader);
	static void calculateSolutionDynamicProgramming(FileReader* fileReader);
	static void calculateSolutionSimulatedAnnealing(FileReader* fileReader);
	static void calculateSolutionGeneticAlgorithm(FileReader* fileReader);
};