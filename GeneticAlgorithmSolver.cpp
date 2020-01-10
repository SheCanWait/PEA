#include "GeneticAlgorithmSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>

GeneticAlgorithmSolver::GeneticAlgorithmSolver(std::vector<Edge> edges, int numberOfNodes)
	: edges(edges), numberOfNodes(numberOfNodes), solutionCost(2147483647) {
	for (int i = 0; i < numberOfNodes; i++)
	{
		nodes.push_back(i);
	}
	convertToMatrix();
}

void GeneticAlgorithmSolver::convertToMatrix() {
	for (int i = 0; i < numberOfNodes; i++) {
		std::vector<int> singleVector;
		for (int j = 0; j < numberOfNodes; j++) {
			singleVector.push_back(edges.at(i * numberOfNodes + j).getCost());
		}
		matrix.push_back(singleVector);
	}
}

void GeneticAlgorithmSolver::printEdges(std::string caseName) {
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

int GeneticAlgorithmSolver::calculateCost(std::vector<int> orderOfCities) {
	int cost = 0;
	for (int i = 0; i < orderOfCities.size() - 1; i++) {
		cost += matrix.at(orderOfCities.at(i + 1)).at(orderOfCities.at(i));
	}

	cost += matrix.at(orderOfCities.at(0)).at(orderOfCities.at(orderOfCities.size() - 1));
	return cost;
}

bool operator== (const SimulatedAnnealingSolution& lhs, const SimulatedAnnealingSolution& rhs) {
	return lhs.cost == rhs.cost;
}

bool operator< (const SimulatedAnnealingSolution& lhs, const SimulatedAnnealingSolution& rhs) {
	return lhs.cost < rhs.cost;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::generateRandomPopulation(int populationSize) {
	int matrixSize = matrix.size();
	std::vector<SimulatedAnnealingSolution> population;
	for (int i = 0 ; i < populationSize ; i++) {
		std::vector<int> path;
		while (path.size() < matrixSize) {
			int element = std::rand() % matrixSize;
			if (std::find(path.begin(), path.end(), element) != path.end()) {
				continue;
			}
			else {
				path.push_back(element);
			}
		}
		SimulatedAnnealingSolution solutionToCompare;
		solutionToCompare.cost = calculateCost(path);
		if (std::find(population.begin(), population.end(), solutionToCompare) != population.end()) {
			i--;
			continue;
		}
		else {
			SimulatedAnnealingSolution individual;
			individual.orderOfCities = path;
			individual.cost = calculateCost(individual.orderOfCities);
			population.push_back(individual);
		}
	}
	return population;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::addCrossingsToPopulation(std::vector<SimulatedAnnealingSolution> population, int numberOfCrossings) {
	for (int i = 0; i < numberOfCrossings; i++) {
		int firstIndividualIndex = std::rand() % population.size();
		int secondIndividualIndex = std::rand() % population.size();
		while (firstIndividualIndex == secondIndividualIndex) {
			secondIndividualIndex = std::rand() % population.size();
		}
		population.push_back(addHalfCrossingFirstChild(population.at(firstIndividualIndex), population.at(secondIndividualIndex)));
		population.push_back(addHalfCrossingSecondChild(population.at(firstIndividualIndex), population.at(secondIndividualIndex)));
	}
	return population;
}

SimulatedAnnealingSolution GeneticAlgorithmSolver::addHalfCrossingFirstChild(SimulatedAnnealingSolution firstIndividual, SimulatedAnnealingSolution secondIndividual) {
	SimulatedAnnealingSolution outputIndividual;
	for (int i = 0; i < firstIndividual.orderOfCities.size() / 2; i++) {
		outputIndividual.orderOfCities.push_back(firstIndividual.orderOfCities.at(i));
	}
	while (outputIndividual.orderOfCities.size() < firstIndividual.orderOfCities.size()) {
		for (int j = 0; j < secondIndividual.orderOfCities.size(); j++) {
			if (std::find(outputIndividual.orderOfCities.begin(), outputIndividual.orderOfCities.end(), secondIndividual.orderOfCities.at(j)) != outputIndividual.orderOfCities.end()) {
				continue;
			}
			else {
				outputIndividual.orderOfCities.push_back(secondIndividual.orderOfCities.at(j));
			}
		}
	}
	outputIndividual.cost = calculateCost(outputIndividual.orderOfCities);
	return outputIndividual;
}

SimulatedAnnealingSolution GeneticAlgorithmSolver::addHalfCrossingSecondChild(SimulatedAnnealingSolution firstIndividual, SimulatedAnnealingSolution secondIndividual) {
	SimulatedAnnealingSolution outputIndividual;
	for (int i = firstIndividual.orderOfCities.size() / 2; i < firstIndividual.orderOfCities.size(); i++) {
		outputIndividual.orderOfCities.push_back(firstIndividual.orderOfCities.at(i));
	}
	while (outputIndividual.orderOfCities.size() < firstIndividual.orderOfCities.size()) {
		for (int j = 0; j < secondIndividual.orderOfCities.size(); j++) {
			if (std::find(outputIndividual.orderOfCities.begin(), outputIndividual.orderOfCities.end(), secondIndividual.orderOfCities.at(j)) != outputIndividual.orderOfCities.end()) {
				continue;
			}
			else {
				outputIndividual.orderOfCities.push_back(secondIndividual.orderOfCities.at(j));
			}
		}
	}
	outputIndividual.cost = calculateCost(outputIndividual.orderOfCities);
	return outputIndividual;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::addCrossingsToPopulationPMX(std::vector<SimulatedAnnealingSolution> population, int numberOfCrossings) {
	for (int i = 0; i < numberOfCrossings; i++) {
		int firstIndividualIndex = std::rand() % population.size();
		int secondIndividualIndex = std::rand() % population.size();
		while (firstIndividualIndex == secondIndividualIndex) {
			secondIndividualIndex = std::rand() % population.size();
		}
		SimulatedAnnealingSolution parent1 = population.at(firstIndividualIndex);
		SimulatedAnnealingSolution parent2 = population.at(secondIndividualIndex);

		std::vector<int> initialCityOrder;
		for (int i = 0; i < population.at(0).orderOfCities.size(); i++) {
			initialCityOrder.push_back(0);
		}
		SimulatedAnnealingSolution child;
		child.orderOfCities = initialCityOrder;
		int firstIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
		int secondIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
		while (firstIndexOfCrossing >= secondIndexOfCrossing) {
			firstIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
			secondIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
		}
		for (int i = firstIndexOfCrossing; i < secondIndexOfCrossing; i++) {
			child.orderOfCities.at(i) = parent1.orderOfCities.at(i);
		}

		std::vector<std::pair<int, int>> pairs;
		for (int i = firstIndexOfCrossing; i < secondIndexOfCrossing; i++) {
			if (!(std::find(child.orderOfCities.begin(), child.orderOfCities.end(), parent2.orderOfCities.at(i)) != child.orderOfCities.end())) {
				//parent2 doesn't contain chosen element on according positions in parent1
				std::pair<int, int> pair;
				pair.first = parent2.orderOfCities.at(i);
				pair.second = child.orderOfCities.at(i);
				pairs.push_back(pair);
			}
		}

		for (int i = 0; i < pairs.size(); i++) {
			std::vector<int>::iterator it = std::find(parent2.orderOfCities.begin(), parent2.orderOfCities.end(), pairs.at(i).second);
			int secondElementInPairInSecondParentIndex = std::distance(parent2.orderOfCities.begin(), it);
			if (child.orderOfCities.at(secondElementInPairInSecondParentIndex) == 0) {
				child.orderOfCities.at(secondElementInPairInSecondParentIndex) = pairs.at(i).first;
			}
			else {
				pairs.at(i).first = pairs.at(i).second;
				std::vector<int>::iterator it2 = std::find(parent1.orderOfCities.begin(), parent1.orderOfCities.end(), pairs.at(i).first);
				int firstElementInPairInFirstParentIndex = std::distance(parent1.orderOfCities.begin(), it2);
				pairs.at(i).second = parent2.orderOfCities.at(firstElementInPairInFirstParentIndex);
				i--;
			}
		}
		population.push_back(child);
		population.push_back(PMXForSecondChild(secondIndividualIndex, firstIndividualIndex, population));
	}
	return population;
}

SimulatedAnnealingSolution GeneticAlgorithmSolver::PMXForSecondChild(int secondParent, int firstParent, std::vector<SimulatedAnnealingSolution> population) {
	int firstIndividualIndex = secondParent;
	int secondIndividualIndex = firstParent;

	SimulatedAnnealingSolution parent1 = population.at(firstIndividualIndex);
	SimulatedAnnealingSolution parent2 = population.at(secondIndividualIndex);

	std::vector<int> initialCityOrder;
	for (int i = 0; i < population.at(0).orderOfCities.size(); i++) {
		initialCityOrder.push_back(0);
	}
	SimulatedAnnealingSolution child;
	child.orderOfCities = initialCityOrder;
	int firstIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
	int secondIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
	while (firstIndexOfCrossing >= secondIndexOfCrossing) {
		firstIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
		secondIndexOfCrossing = std::rand() % parent1.orderOfCities.size();
	}
	for (int i = firstIndexOfCrossing; i < secondIndexOfCrossing; i++) {
		child.orderOfCities.at(i) = parent1.orderOfCities.at(i);
	}

	std::vector<std::pair<int, int>> pairs;
	for (int i = firstIndexOfCrossing; i < secondIndexOfCrossing; i++) {
		if (!(std::find(child.orderOfCities.begin(), child.orderOfCities.end(), parent2.orderOfCities.at(i)) != child.orderOfCities.end())) {
			//parent2 doesn't contain chosen element on according positions in parent1
			std::pair<int, int> pair;
			pair.first = parent2.orderOfCities.at(i);
			pair.second = child.orderOfCities.at(i);
			pairs.push_back(pair);
		}
	}

	for (int i = 0; i < pairs.size(); i++) {
		std::vector<int>::iterator it = std::find(parent2.orderOfCities.begin(), parent2.orderOfCities.end(), pairs.at(i).second);
		int secondElementInPairInSecondParentIndex = std::distance(parent2.orderOfCities.begin(), it);
		if (child.orderOfCities.at(secondElementInPairInSecondParentIndex) == 0) {
			child.orderOfCities.at(secondElementInPairInSecondParentIndex) = pairs.at(i).first;
		}
		else {
			pairs.at(i).first = pairs.at(i).second;
			std::vector<int>::iterator it2 = std::find(parent1.orderOfCities.begin(), parent1.orderOfCities.end(), pairs.at(i).first);
			int firstElementInPairInFirstParentIndex = std::distance(parent1.orderOfCities.begin(), it2);
			pairs.at(i).second = parent2.orderOfCities.at(firstElementInPairInFirstParentIndex);
			i--;
		}
	}
	return child;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::createMatingPoolByRanking(std::vector<SimulatedAnnealingSolution> population, int populationSize) {
	std::sort(population.begin(), population.end());
	std::vector<SimulatedAnnealingSolution> outputPopulation;
	for (int i = 0; i < populationSize; i++) {
		outputPopulation.push_back(population.at(i));
	}
	std::sort(outputPopulation.begin(), outputPopulation.end());
	return outputPopulation;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::createMatingPoolByRoulette(std::vector<SimulatedAnnealingSolution> population, int populationSize) {
	SimulatedAnnealingSolution worstSolution = population.at(0);
	for (int i = 0; i < population.size(); i++) {
		if (population.at(i).cost > worstSolution.cost) {
			worstSolution = population.at(i);
		}
	}

	int sumOfFitnessOfPopulation = 0;
	for (int i = 0; i < population.size(); i++) {
		sumOfFitnessOfPopulation += population.at(i).cost;
	}

	std::vector<SimulatedAnnealingSolution> vectorOfProbability;
	for (int i = 0; i < population.size(); i++) {
		for (int j = 0; j < worstSolution.cost - population.at(i).cost + 1; j++) {
			vectorOfProbability.push_back(population.at(i));
		}
	}
	std::sort(vectorOfProbability.begin(), vectorOfProbability.end());
	std::vector<SimulatedAnnealingSolution> outputPopulation;
	for (int i = 0; i < populationSize; i++) {
		int randomVectorIndex = std::rand() % vectorOfProbability.size();
		outputPopulation.push_back(vectorOfProbability.at(randomVectorIndex));
		vectorOfProbability.erase(std::remove_if(vectorOfProbability.begin(), vectorOfProbability.end(),
			[outputPopulation](SimulatedAnnealingSolution const& x) { return x.cost == outputPopulation.at(outputPopulation.size() - 1).cost 
			&& x.orderOfCities == outputPopulation.at(outputPopulation.size() - 1).orderOfCities; 
			}), vectorOfProbability.end());
	}
	std::sort(outputPopulation.begin(), outputPopulation.end());
	return outputPopulation;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::mutatePopulationBySwap(std::vector<SimulatedAnnealingSolution> population) {
	int individualIndex = std::rand() % population.size();
	std::vector<int> citiesInOrder = population.at(individualIndex).orderOfCities;
	population.erase(population.begin() + individualIndex);

	int firstCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	int secondCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	while (firstCityToSwap == secondCityToSwap) {
		secondCityToSwap = std::rand() % (citiesInOrder.size() - 1);
	}
	int firstCity = citiesInOrder.at(firstCityToSwap);
	int secondCity = citiesInOrder.at(secondCityToSwap);
	citiesInOrder.at(firstCityToSwap) = secondCity;
	citiesInOrder.at(secondCityToSwap) = firstCity;

	SimulatedAnnealingSolution individual;
	individual.orderOfCities = citiesInOrder;
	individual.cost = calculateCost(individual.orderOfCities);
	population.push_back(individual);
	return population;
}

std::vector<SimulatedAnnealingSolution> GeneticAlgorithmSolver::mutatePopulationByInversion(std::vector<SimulatedAnnealingSolution> population) {
	int individualIndex = std::rand() % population.size();
	std::vector<int> citiesInOrder = population.at(individualIndex).orderOfCities;
	population.erase(population.begin() + individualIndex);

	int inversionBeginningIndex = std::rand() % (citiesInOrder.size() - 1);
	int inversionEndingIndex = std::rand() % (citiesInOrder.size() - 1);
	while (inversionBeginningIndex >= inversionEndingIndex) {
		inversionBeginningIndex = std::rand() % (citiesInOrder.size() - 1);
		inversionEndingIndex = std::rand() % (citiesInOrder.size() - 1);
	}
	std::reverse(citiesInOrder.begin() + inversionBeginningIndex, citiesInOrder.begin() + inversionEndingIndex);

	SimulatedAnnealingSolution individual;
	individual.orderOfCities = citiesInOrder;
	individual.cost = calculateCost(individual.orderOfCities);
	population.push_back(individual);
	return population;
}

float GeneticAlgorithmSolver::getRandomNumberFrom0To1() {
	float number = std::rand();
	number = number / (float)RAND_MAX;
	return number;
}

SimulatedAnnealingSolution GeneticAlgorithmSolver::solve(int numberOfRepetitions, int populationSize, int numberOfCrossings, int eliteSize, float mutationRate, int crossingRate) {
	SimulatedAnnealingSolution finalSolution;
	finalSolution.cost = INT_MAX;

	std::vector<SimulatedAnnealingSolution> initialPopulation = generateRandomPopulation(populationSize);
	std::vector<SimulatedAnnealingSolution> matingPoolPopulation = initialPopulation;

	for (int i = 0; i < numberOfRepetitions; i++) {
		SimulatedAnnealingSolution currentPopulationSolution;
		currentPopulationSolution = matingPoolPopulation.at(0);
		if (currentPopulationSolution.cost < finalSolution.cost) {
			finalSolution = currentPopulationSolution;
		}

		std::vector<SimulatedAnnealingSolution> populationWithCrossings = matingPoolPopulation;
		if (getRandomNumberFrom0To1() < crossingRate) {
			populationWithCrossings = addCrossingsToPopulationPMX(matingPoolPopulation, numberOfCrossings);
		}

		std::vector<SimulatedAnnealingSolution> mutatedPopulation = populationWithCrossings;
		if (getRandomNumberFrom0To1() < mutationRate) {
			mutatedPopulation = mutatePopulationBySwap(populationWithCrossings);
		}

		matingPoolPopulation = createMatingPoolByRanking(mutatedPopulation, eliteSize);
	}
	return finalSolution;
}
