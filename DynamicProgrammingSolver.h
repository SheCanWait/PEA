#pragma once
#include <vector>
#include <string>
#include "Edge.h"
#include "Solution.h"
#include "SearchWrapper.h"
#include <unordered_map>
#include <map>

class DynamicProgrammingSolver
{
public:
	DynamicProgrammingSolver(std::vector<Edge> edges, int numberOfNodes);
	void printEdges(std::string caseName);
	void printSolution(Solution solution);
	void convertToMatrix();
	std::vector<int> optimalOrderOfCities;
	std::vector<int> tempOrder;
	Solution solveRecursive(int mask, int startCity);
	std::map<SearchWrapper, Solution> solutionsMap;
private:
	std::vector<std::vector<int>> matrix;
	std::vector<Edge> edges;
	std::vector<int> nodes;
	int numberOfNodes;
	int solutionCost;
	std::string caseName;
};

