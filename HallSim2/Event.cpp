#include "stdafx.h"
#include "Event.h"
#include "Cell.h"
#include <stdlib.h>

namespace std {
 
Event::Event()
{
	//Cell currentCell;
	//int timePoint = rand() % 6 + 5;
	currentCell = new Cell();
	timePoint = rand() % 6 + 5;
}


//Event::Event(Cell * x, int t){
//	currentCell = *x;
//	timePoint = t;
//}
//
//Event::Event(CancerCell * x, int t){
//	currentCell = *x;
//	timePoint = t;
//}
//
//void Event::setCell(Cell & newCell){
//	currentCell = newCell;
//}
//
//Cell * Event::getCell(){
//	Cell *temp = &currentCell;
//	return temp;
//	//return currentCell; 
//}

//Second try.. pass by reference... pg 38 
Event::Event(Cell & x, int t){
	currentCell = &x;
	timePoint = t + rand() % 6 + 5;
}

/**
Create the new event with the cell passed it
Schedule it's mitosis for some point in time
The time point is a bit longer if it is on the diagonal
*/
Event::Event(Cell & x, int t, int dir){
	if (dir == 4 || dir == 5 || dir==6 || dir==7){
		timePoint = t + rand() % 7 + 7;
		//Range from 7 to 14 time steps in future
		//(sqrt(2)) of 5 to 11 time steps to account for 
		//diagonal taking more time
	}
	else {
		timePoint = t + rand() % 6 + 5;
		//Range from 5 to 11 time steps in future
	}
	currentCell = &x;
}

//Event::Event(CancerCell & x, int t){
//	currentCell = &x;
//	timePoint = t + rand() % 6 + 5;
//}

void Event::setCell(Cell & newCell){
	currentCell = &newCell;
}

Cell * Event::getCell(){
	return currentCell;
	//return currentCell; 
}

void Event::setTime(int newTime){
	timePoint = newTime;
}

Event::~Event(){
	//TO-DO make deconstructor
}

}