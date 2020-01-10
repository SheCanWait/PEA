#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "Solver.h"
#include "Solution.h"

Solver::Solver(std::vector<Edge> edges, int numberOfNodes)
	: edges(edges), numberOfNodes(numberOfNodes), solutionCost(2147483647) {
	for (int i = 0; i < numberOfNodes; i++)
	{
		nodes.push_back(i);
	}
}

void Solver::printEdges(std::string caseName) {
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

Solution Solver::getSolutionForInstance() {
	Solution* solution = new Solution();
	solution->citiesInOrder = *new std::vector<int>;
	solution->cost = INT_MAX;
	std::vector<int> citiesInOrder = createOrderOfCities();
	do {
		int currentCost = calculateCost(citiesInOrder);
		if (currentCost < solution->cost) {
			std::cout << currentCost << std::endl;
			solution->citiesInOrder = citiesInOrder;
			solution->cost = currentCost;
		}
	} while (std::next_permutation(&(citiesInOrder.at(0)), &(citiesInOrder.at(numberOfNodes - 1))));
	return *solution;
}

std::vector<int> Solver::createOrderOfCities() {
	std::vector<int> citiesInOrder;
	for (int i = 0; i < numberOfNodes; i++) {
		citiesInOrder.push_back(i);
	}
	return citiesInOrder;
}

int Solver::calculateCost(std::vector<int> citiesInOrder) {
	int cost = 0;
	for (int i = 0 ; i < citiesInOrder.size() - 1 ; i++) {
		int source = citiesInOrder[i];
		int destination = citiesInOrder[i + 1];
		cost += edges.at(source * numberOfNodes + destination).getCost();
	}
	cost += edges.at(citiesInOrder.at(citiesInOrder.size() - 1) * numberOfNodes + citiesInOrder.at(0)).getCost();
	return cost;
}

void Solver::printSolution(Solution solution) {
	for (int i = 0; i < numberOfNodes - 1; i++) {
		std::cout << solution.citiesInOrder[i] << "->";
	}
	std::cout << solution.citiesInOrder[numberOfNodes - 1];
	std::cout << "\nSolution cost: " << solution.cost << "\n";
}