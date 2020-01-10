#include "FileReader.h"

FileReader::FileReader() {}

FileReader::FileReader(std::string fileName)
	: inputFileName(fileName) {
	readInput();
}

std::vector<Edge> FileReader::getEdges() const {
	return edges;
}

int FileReader::getNumberOfNodes() const {
	return numberOfEdgesInLine;
}

std::string FileReader::getCaseName() const {
	return caseName;
}

bool FileReader::readInput() {
	inputFile.open(inputFileName);
	if (!inputFile)
	{
		return false;
	}
	std::getline(inputFile, caseName);
	inputFile >> numberOfEdgesInLine;
	for (int i = 0; i < numberOfEdgesInLine; i++)
	{
		for (int j = 0; j < numberOfEdgesInLine; j++)
		{
			int cost;
			inputFile >> cost;
			Edge edge(i, j, cost);
			edges.push_back(edge);
		}
	}
	return true;
}
