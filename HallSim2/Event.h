#ifndef EVENT_H_
#define EVENT_H_

#include "Cell.h"	
 
namespace std {

//I'm trying to create an event class that is just the cell, and how many timesteps to add to the counter in the event
	//the next item in the event queue is not one time step ahead.. I want to be able to add the right number of steps to the queue
	//This thing should take a Cell or a CancerCell... since CancerCell is a child of Cell, shouldn't that be easy to do?
 
//class Event {
//public:
//	Event();
//	Event(Cell * x, int t);
//	Event(CancerCell * x, int t);
//	int getTime() { return timePoint; }
//	void setTime(int newTime);
//	Cell * getCell(); //TO-DO is this being passed by value or reference...??
//	void setCell(Cell & newCell);
//	virtual ~Event();
//private:
//	Cell currentCell;
//	int timePoint;
//};

	//Second try... just using the actual objects
	//I think I want the events to just hold points to Cells, then using polymorphism and virtual members I can send in cancer cells too ...
	//Trying to get it working with regular cells first 
class Event {
public:
	Event();
	Event(Cell & x, int t); //pass by reference
	//Event(CancerCell & x, int t); //pass by reference, pg 38
	int getTime() { return timePoint; }
	void setTime(int newTime); //look up what an ampersand means in a function prototype
	Cell * getCell(); //TO-DO is this being passed by value or reference...?? By reference I think
	void setCell(Cell & newCell);
	virtual ~Event();
private:
	Cell * currentCell; //Class is made up of a pointer to a cell and a time
	int timePoint;
};
 
 
} /* namespace std */
#endif /* EVENT_H_ */