#pragma once
#include <vector>
#include "FileReader.h"
#include "Edge.h"
#include "Solution.h"
#include <string>

class Solver
{
public:
	Solver(std::vector<Edge> edges, int numberOfNodes);
	void printEdges(std::string caseName);
	Solution getSolutionForInstance();
	void printSolution(Solution solution);
private:
	std::vector<int> createOrderOfCities();
	int calculateCost(std::vector<int> citiesInOrder);
	std::vector<Edge> edges;
	std::vector<int> nodes;
	std::vector<int> solutions;
	int numberOfNodes;
	int solutionCost;
	std::string caseName;
};

