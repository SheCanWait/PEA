#pragma once
#include <string>
class Edge
{
public:
	Edge(int origin, int destination, int cost)
		: originID(origin), destinationID(destination), cost(cost) {};

	int getCost() const { return cost; };
	int getDestinationID() const { return destinationID; };
	int getOriginID() const { return originID; };
private:
	int originID;
	int destinationID;
	int cost;
};

