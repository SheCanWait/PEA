#include <iostream>
#include <string>
#include "PEA_0.h"
#include "FileReader.h"
#include "Solver.h"
#include "DynamicProgrammingSolver.h"
#include <chrono>
#include "SimulatedAnnealingSolver.h"
#include "SimulatedAnnealingSolution.h"
#include "GeneticAlgorithmSolver.h"

int main() {
	srand(time(NULL));
	std::string fileNames[9] = {"data10.txt", "data11.txt", "data12.txt", "data13.txt", "data14.txt", "data15.txt", "data16.txt", "data17.txt", "data18.txt"};
	std::string fineNamesATSP[18] = { "data17.txt", "data34.txt", "data36.txt", "data39.txt", "data43.txt", "data45.txt", "data48.txt",  "data53.txt",  "data56.txt",  "data65.txt",  "data70.txt",  "data71.txt",  "data100.txt",  "data171.txt",  "data323.txt",  "data358.txt",  "data403.txt",  "data443.txt" };


	FileReader* fileReader = new FileReader();
	PEA_0 pea0;

	for (std::string fileName : fileNames) {
		fileReader = new FileReader(fileName);
		auto time1 = std::chrono::high_resolution_clock::now();
		//pea0.calculateSolution(fileReader);
		pea0.calculateSolutionGeneticAlgorithm(fileReader);
		auto time2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = time2 - time1;
		std::cout << "time: " << diff.count() << "\n";
	}
	for (std::string fileName : fineNamesATSP) {
		fileReader = new FileReader(fileName);
		auto time1 = std::chrono::high_resolution_clock::now();
		//pea0.calculateSolution(fileReader);
		pea0.calculateSolutionGeneticAlgorithm(fileReader);
		auto time2 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> diff = time2 - time1;
		std::cout << "time: " << diff.count() << "\n";
	}
	return 0;
}

void PEA_0::calculateSolution(FileReader* fileReader) {
	Solver solver((*fileReader).getEdges(), (*fileReader).getNumberOfNodes());
	solver.printEdges((*fileReader).getCaseName());
	Solution solution = solver.getSolutionForInstance();
	solver.printSolution(solution);
}

void PEA_0::calculateSolutionDynamicProgramming(FileReader* fileReader) {
	DynamicProgrammingSolver solver((*fileReader).getEdges(), (*fileReader).getNumberOfNodes());
	solver.printEdges((*fileReader).getCaseName());
	std::vector<int> vector;
	Solution solution = solver.solveRecursive(1, 0);
	std::cout << "optimal cost = " << solution.cost << std::endl;
	for (int i = 0; i < solution.citiesInOrder.size() - 1; i++) {
		std::cout << solution.citiesInOrder.at(i) << "->";
	}
	std::cout << solution.citiesInOrder.at(solution.citiesInOrder.size() - 1) << std::endl;
}

void PEA_0::calculateSolutionSimulatedAnnealing(FileReader* fileReader) {
	SimulatedAnnealingSolver solver((*fileReader).getEdges(), (*fileReader).getNumberOfNodes());

	float initialTemperature = 10000;
	float endingTemperature = 10;
	float const alpha = 0.995;
	int repetitionsForOneTemperature = 50;

	SimulatedAnnealingSolution solution = solver.solve(initialTemperature, endingTemperature, alpha, repetitionsForOneTemperature);
	std::cout << fileReader->getCaseName() << std::endl;
	std::cout << "optimal cost = " << solution.cost << std::endl;
	for (int i = 0; i < solution.orderOfCities.size() - 1; i++) {
		std::cout << solution.orderOfCities.at(i) << "->";
	}
	std::cout << solution.orderOfCities.at(solution.orderOfCities.size() - 1) << std::endl;
}

void PEA_0::calculateSolutionGeneticAlgorithm(FileReader* fileReader) {
	GeneticAlgorithmSolver solver((*fileReader).getEdges(), (*fileReader).getNumberOfNodes());

	int populationSize = 300;
	int numberOfCrossings = 800;
	int numberOfRepetitions = 1000;
	int eliteSize = 300; //has to be <= populationSize
	float mutationRate = 0.05;
	float crossingRate = 0.8;

	SimulatedAnnealingSolution solution = solver.solve(numberOfRepetitions, populationSize, numberOfCrossings, eliteSize, mutationRate, crossingRate);
	std::cout << fileReader->getCaseName() << std::endl;
	std::cout << "optimal cost = " << solution.cost << ", ";
}