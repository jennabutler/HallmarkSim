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