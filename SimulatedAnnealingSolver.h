#pragma once
#include <vector>
#include <string>
#include "Edge.h"
#include "SimulatedAnnealingSolution.h"

class SimulatedAnnealingSolver
{
public:

	std::vector<std::vector<int>> matrix;
	std::string caseName;
	int numberOfNodes;
	int solutionCost;
	std::vector<Edge> edges;
	std::vector<int> nodes;

	SimulatedAnnealingSolver(std::vector<Edge> edges, int numberOfNodes);
	void convertToMatrix();
	void printEdges(std::string caseName);
	int calculateCost(std::vector<int> orderOfCities);
	static float getrandomNumberFrom0To1();

	float temperature;
	float endingTemperature;
	
	SimulatedAnnealingSolution solve(float initialTemperature, float endingTemperature, float alpha, int repetitionsForOneTemperature);
	std::vector<int> getNeighbourBySwap(std::vector<int> citiesInOrder);
	std::vector<int> getNeighbourByReverse(std::vector<int> citiesInOrder);
	float getTemperature(int currentTemperature, float alpha);
	float getProbability(float e, float en);
	float getCriterion(SimulatedAnnealingSolution solution);
};

