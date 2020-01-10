#pragma once
#include <fstream>
#include <string>
#include <vector>
#include "Edge.h"

class FileReader
{
public:
	FileReader();
	FileReader(std::string fileName);
	std::vector<Edge> getEdges() const;
	int getNumberOfNodes() const;
	std::string getCaseName() const;
private:
	std::ifstream inputFile;
	std::string inputFileName;
	std::string caseName;
	int numberOfEdgesInLine;
	std::vector<Edge> edges;

	bool readInput();
};

