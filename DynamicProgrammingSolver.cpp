#include "DynamicProgrammingSolver.h"
#include <iostream>
#include <iomanip>
#include <cmath>

DynamicProgrammingSolver::DynamicProgrammingSolver(std::vector<Edge> edges, int numberOfNodes)
	: edges(edges), numberOfNodes(numberOfNodes), solutionCost(2147483647) {
		for (int i = 0; i < numberOfNodes; i++)
		{
			nodes.push_back(i);
		}
		convertToMatrix();
}

void DynamicProgrammingSolver::convertToMatrix() {
	for (int i = 0; i < numberOfNodes; i++) {
		std::vector<int> singleVector;
		for (int j = 0; j < numberOfNodes; j++) {
			singleVector.push_back(edges.at(i * numberOfNodes + j).getCost());
		}
		matrix.push_back(singleVector);
	}
}

void DynamicProgrammingSolver::printEdges(std::string caseName) {
	int i = 0;
	std::cout << caseName << "\n";
	for (Edge edge : edges)
	{
		std::cout << std::setw(4) << std::setfill(' ') << edge.getCost() << " ";
		i++;
		if (i == numberOfNodes)
		{
			i = 0;
			std::cout << "\n";
		}
	}
}

void DynamicProgrammingSolver::printSolution(Solution solution) {
	for (int i = 0 ; i < numberOfNodes - 1 ; i++) {
		std::cout << solution.citiesInOrder[i] << "->";
	}
	std::cout << solution.citiesInOrder[numberOfNodes - 1];
	std::cout << "\nSolution cost: " << solution.cost << "\n";
}



Solution DynamicProgrammingSolver::solveRecursive(int mask, int startCity) {
	Solution solution1;
	if (mask == (1 << matrix.size()) - 1) {
		solution1.citiesInOrder.push_back(startCity);
		solution1.cost = matrix[startCity][0];
		return solution1;
	}
	Solution solution2;
	solution2.cost = INT_MAX;

	for (int city = 0; city < matrix.size(); city++) {
		if ((mask & (1 << city)) == 0) {
			SearchWrapper searchWrapper;
			searchWrapper.startCity = city;
			searchWrapper.mask = mask | (1 << city);
			Solution solutionFromMap = solutionsMap[searchWrapper];
			Solution someSolution;
			if (solutionFromMap.cost == 0) {
				someSolution = solveRecursive(mask | (1 << city), city);
				SearchWrapper searchWrapperToInsert;
				searchWrapperToInsert.startCity = city;
				searchWrapperToInsert.mask = mask | (1 << city);
				solutionsMap[searchWrapperToInsert] = someSolution;
			}
			else {
				someSolution = solutionFromMap;
			}
			Solution finalSolution;
			finalSolution.cost = matrix[startCity][city] + someSolution.cost;
			finalSolution.citiesInOrder.push_back(startCity);
			finalSolution.citiesInOrder.insert(finalSolution.citiesInOrder.end(),
				someSolution.citiesInOrder.begin(), someSolution.citiesInOrder.end());
			if (solution2.cost > finalSolution.cost) {
				solution2 = finalSolution;
			}
		}
	}
	return solution2;
}

bool operator< (const SearchWrapper& lhs, const SearchWrapper& rhs) {
	if (lhs.mask < rhs.mask && lhs.startCity < rhs.startCity) {
		return true;
	}
	if (lhs.mask == rhs.mask && lhs.startCity == rhs.startCity) {
		return false;
	}
	if (lhs.mask > rhs.mask && lhs.startCity > rhs.startCity) {
		return false;
	}
	if (lhs.mask > rhs.mask && lhs.startCity < rhs.startCity) {
		return true;
	}
	if (lhs.mask < rhs.mask && lhs.startCity > rhs.startCity) {
		return false;
	}
	if (lhs.mask == rhs.mask && lhs.startCity < rhs.startCity) {
		return true;
	}
	if (lhs.mask == rhs.mask && lhs.startCity > rhs.startCity) {
		return false;
	}
	if (lhs.mask > rhs.mask && lhs.startCity == rhs.startCity) {
		return false;
	}
	if (lhs.mask < rhs.mask && lhs.startCity == rhs.startCity) {
		return true;
	}
}
