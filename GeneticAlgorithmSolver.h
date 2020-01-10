#pragma once
#include <vector>
#include <string>
#include "Edge.h"
#include "SimulatedAnnealingSolution.h"

class GeneticAlgorithmSolver
{
public:
	std::vector<std::vector<int>> matrix;
	std::string caseName;
	int numberOfNodes;
	int solutionCost;
	std::vector<Edge> edges;
	std::vector<int> nodes;

	GeneticAlgorithmSolver(std::vector<Edge> edges, int numberOfNodes);
	void convertToMatrix();
	void printEdges(std::string caseName);
	int calculateCost(std::vector<int> orderOfCities);

	SimulatedAnnealingSolution solve(int numberOfRepetitions, int populationSize, int numberOfCrossings, int eliteSize, float mutationRate, int crossingRate);
	std::vector<SimulatedAnnealingSolution> generateRandomPopulation(int populationSize);
	std::vector<SimulatedAnnealingSolution> addCrossingsToPopulation(std::vector<SimulatedAnnealingSolution> population, int numberOfCrossings);
	std::vector<SimulatedAnnealingSolution> addCrossingsToPopulationPMX(std::vector<SimulatedAnnealingSolution> population, int numberOfCrossings);
	SimulatedAnnealingSolution PMXForSecondChild(int secondIndividualIndex, int firstIndividualIndex, std::vector<SimulatedAnnealingSolution> population);
	SimulatedAnnealingSolution addHalfCrossingFirstChild(SimulatedAnnealingSolution firstIndividual, SimulatedAnnealingSolution secondIndividual);
	SimulatedAnnealingSolution addHalfCrossingSecondChild(SimulatedAnnealingSolution firstIndividual, SimulatedAnnealingSolution secondIndividual);
	std::vector<SimulatedAnnealingSolution> mutatePopulationBySwap(std::vector<SimulatedAnnealingSolution> population);
	std::vector<SimulatedAnnealingSolution> mutatePopulationByInversion(std::vector<SimulatedAnnealingSolution> population);
	std::vector<SimulatedAnnealingSolution> createMatingPoolByRanking(std::vector<SimulatedAnnealingSolution> population, int populationSize);
	std::vector<SimulatedAnnealingSolution> createMatingPoolByRoulette(std::vector<SimulatedAnnealingSolution> population, int populationSize);
	float getRandomNumberFrom0To1();
};

