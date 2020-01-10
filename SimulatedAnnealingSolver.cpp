#include "SimulatedAnnealingSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

const float EULER = 2.71828182845904523536;


SimulatedAnnealingSolver::SimulatedAnnealingSolver(std::vector<Edge> edges, int numberOfNodes)
	: edges(edges), numberOfNodes(numberOfNodes), solutionCost(2147483647) {
	for (int i = 0; i < numberOfNodes; i++)
	{
		nodes.push_back(i);
	}
	convertToMatrix();
}

void SimulatedAnnealingSolver::convertToMatrix() {
	for (int i = 0; i < numberOfNodes; i++) {
		std::vector<int> singleVector;
		for (int j = 0; j < numberOfNodes; j++) {
			singleVector.push_back(edges.at(i * numberOfNodes + j).getCost());
		}
		matrix.push_back(singleVector);
	}
}

void SimulatedAnnealingSolver::printEdges(std::string caseName) {
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

int SimulatedAnnealingSolver::calculateCost(std::vector<int> orderOfCities) {
	int cost = 0;
	for (int i = 0; i < orderOfCities.size() - 1; i++) {
		cost += matrix.at(orderOfCities.at(i + 1)).at(orderOfCities.at(i));
	}

	cost += matrix.at(orderOfCities.at(0)).at(orderOfCities.at(orderOfCities.size() - 1));
	return cost;
}

std::vector<int> SimulatedAnnealingSolver::getNeighbourBySwap(std::vector<int> citiesInOrder) {
	int firstCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	int secondCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	while (firstCityToSwap == secondCityToSwap) {
		secondCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	}
	int firstCity = citiesInOrder.at(firstCityToSwap);
	int secondCity = citiesInOrder.at(secondCityToSwap);
	citiesInOrder.at(firstCityToSwap) = secondCity;
	citiesInOrder.at(secondCityToSwap) = firstCity;

	return citiesInOrder;
}

std::vector<int> SimulatedAnnealingSolver::getNeighbourByReverse(std::vector<int> citiesInOrder) {
	int startReverseCity = std::rand() % (citiesInOrder.size() - 1);
	int endReverseCity = std::rand() % (citiesInOrder.size() - 1);
	while (startReverseCity >= endReverseCity) {
		startReverseCity = std::rand() % (citiesInOrder.size() - 1);
		endReverseCity = std::rand() % (citiesInOrder.size() - 1);
	}

	std::vector<int> citiesToRevert;
	for(int i = startReverseCity ; i <= endReverseCity ; i++) {
		citiesToRevert.push_back(citiesInOrder.at(i));
	}

	std::reverse(citiesToRevert.begin(), citiesToRevert.end());
	for (int i = startReverseCity; i <= endReverseCity; i++) {
		citiesInOrder.at(i) = citiesToRevert.at(i - startReverseCity);
	}

	return citiesInOrder;
}

float SimulatedAnnealingSolver::getTemperature(int currentTemperature, float alpha) {
	return currentTemperature * alpha;
}

float SimulatedAnnealingSolver::getProbability(float e, float en) {
	if (e - en < 0) {
		float probability = 1 / (1 + pow(EULER, ((e - en) / temperature)));
		return probability;
	}
	else {
		return 1;
	}
}

float SimulatedAnnealingSolver::getrandomNumberFrom0To1() {
	float number = std::rand();
	number = number / (float)RAND_MAX;
	return number;
}

float SimulatedAnnealingSolver::getCriterion(SimulatedAnnealingSolution solution) {
	return solution.cost;
}

SimulatedAnnealingSolution SimulatedAnnealingSolver::solve(float initialTemperature, float endTemperature, float alpha, int repetitionsForOneTemperature) {
	temperature = initialTemperature;
	endingTemperature = endTemperature;
	
	std::vector<int> randomOrder;
	for (int i = 0; i < matrix.size(); i++) {
		randomOrder.push_back(i);
	}
	int bestGlobalCost = INT_MAX;

	SimulatedAnnealingSolution solutionA;
	solutionA.cost = this->calculateCost(randomOrder);
	solutionA.orderOfCities = randomOrder;

	while (temperature > endingTemperature) {
		int repetitionNumber = 0;
		int bestLocalCost = INT_MAX;
		while (repetitionNumber < repetitionsForOneTemperature) {
			SimulatedAnnealingSolution solutionB;
			solutionB.orderOfCities = getNeighbourByReverse(solutionA.orderOfCities);
			solutionB.cost = calculateCost(solutionB.orderOfCities);
			if (getCriterion(solutionA) < getCriterion(solutionB)) {
				solutionA = solutionB;
			}
			else if (getrandomNumberFrom0To1() < getProbability(getCriterion(solutionB), getCriterion(solutionA))) {
				solutionA = solutionB;
			}
			if (solutionA.cost < bestGlobalCost) {
				bestGlobalCost = solutionA.cost;
			}

			repetitionNumber++;

			if (solutionB.cost < bestLocalCost) {
				bestLocalCost = solutionB.cost;
			}
		}
		temperature = getTemperature(temperature, alpha);
	}
	solutionA.cost = bestGlobalCost;
	return solutionA;
}