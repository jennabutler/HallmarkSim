// HallSim2.cpp : Defines the entry point for the console application.
//
//Testing commits
#include "stdafx.h"
#include "Cell.h"
#include "Event.h"
#include "Constants.h"
#include "Global.h"
#include "BinaryFluid.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <sstream>
#include <sys/types.h>
#include <sys/timeb.h>
#include <ctime>
#include <direct.h>

using namespace std;


/*************************************************************************
**************************************************************************
************************ TO BE IMPLEMENTED TO-DO**************************
//Next steps:
//The boundaires aren't taken care of.. need to have null pointers or something in the
neighbour lsits along the edges so that we can't move off the grid and go into -1 space

No advantage for hallmark 8 (avoids Immunity) because we are not modelling an immune response

TO-DO: If I re run..
2) Change the space boundary to 95%
3) Increase number of simulations
***************************************************************************
***************************************************************************/

enum State {
	ALIVE, DEAD, QUIS, APOP, NEC, IM_DEAD, AGG1, AGG2, AGG3, GLY
};

/**
* Overload the () operator to work for a priority queue
* We want a reverse queue.. the smallest timepoint is at the top
**/

class CompareEvent {
public:
	bool operator()(Event& t1, Event& t2)
	{
		if (t1.getTime() > t2.getTime()) return true;
		return false;
	}
};

/**
This structure holds the combinations of hallmarks for knocking two out
Each has a unique string code
*/
struct Combo {
	string code;
	bool h1;
	bool h2;
	bool h3;
	bool h4;
	bool h5;
	bool h6;
	//bool h7;
	//bool h8;
	bool e2;
	bool em1;
	bool em2;
};

enum Direction {
	N, S, E, W, NE, SE, SW, NW
};

void WriteParams(int howEnd){
	stringstream fileName;
	fileName << "..\\Output\\Working\\Parameters.txt";
	string fn = fileName.str();
	stringstream istream;

	istream << "Ended because: " << howEnd << " (0 is time ran out or out of events, 1 is 95% cancer) \n";
	istream << "Seed: " << RAND_SEED << "\n";
	istream << "GRID_SIDE: " << GRID_SIDE << "\n";
	istream << "Blood radius: " << BLOOD_RADIUS << "\n";
	istream << "Blood radius 2: " << BLOOD_RADIUS_2 << "\n";
	istream << "Blood radius 3: " << BLOOD_RADIUS_3 << "\n";
	istream << "Blood radius 4: " << BLOOD_RADIUS_4 << "\n";
	istream << "Blood H: " << BLOOD_H << "\n";
	istream << "Blood K: " << BLOOD_K << "\n";
	istream << "Growth radius: " << GROWTH_RADIUS << "\n";
	istream << "Growth H: " << GROWTH_H << "\n";
	istream << "Growth K: " << GROWTH_K << "\n";
	istream << "Initial mutation rate: " << MUTATION_RATE << "\n";
	istream << "Oxygen consumption of healthy cells: " << OXYGEN_CON_HEALTHY << "\n";
	istream << "Oxygen consumption of quis cells: " << OXYGEN_CON_QUIS << "\n";
	istream << "Oxygen consumption of dereg cells: " << OXYGEN_DEREG << "\n";
	istream << "Oxygen consumption of agg1 cells: " << OXYGEN_CON_CANCER_AGG1 << "\n";
	istream << "Oxygen consumption of agg2 cells: " << OXYGEN_CON_CANCER_AGG2 << "\n";
	istream << "Oxygen consumption of agg3 cells: " << OXYGEN_CON_CANCER_AGG3 << "\n";


	//Write the cells to a file
	istream << "\n";
	ofstream myfile;
	myfile.open(fn);
	myfile << istream.str();
	myfile.close();
}

/**
* Code to write the current state out to a file
*/
void WriteToFile(std::list<Cell *>::const_iterator iterator, int time, std::list<Cell *>::const_iterator iteratorEnd, string paramSet, int it, int currentTime){
	stringstream fileName;

	fileName << "..\\Output\\Working" << "\\Cells_" << paramSet <<"_01_2_" << RUN << "_" << time << "_it_" << it <<".txt";
	string fn = fileName.str();
	stringstream istream;

	//istream << "i, j, state, mut, sg, igi, aa, it, gu \n";
	//Iterate over the cells and get list to draw
	for (; iterator != iteratorEnd; ++iterator) {
		//(*iterator)->print();
//Run and test is quis to see if taht works, then change python to colour differently if that bit is on
		istream << (*iterator)->geti() << " " << (*iterator)->getj() << " " << (*iterator)->getState() << " " << (*iterator)->isQuis() << " " << (*iterator)->isNec() << " " << (*iterator)->isApop() << " ";
		istream << (*iterator)->isAgg1() << " " << (*iterator)->isAgg2() << " " << (*iterator)->isAgg3() << " " << (*iterator)->isAlive() << " " << (*iterator)->isMutated() << " ";
		istream << (*iterator)->selfGrows() << " " << (*iterator)->ignoresGrowthInhibition() << " " << (*iterator)->avoidsApoptosis() << " ";
		istream << (*iterator)->ignoresTelomere() << " " << (*iterator)->angiogenesis() << " "  << (*iterator)->genomeUnstable() << " " << (*iterator)->avoidsImmunity() << " " << (*iterator)->glycoPheno()  <<"\n";
	}

	//Write the cells to a file
	istream << "\n";
	ofstream myfile;
	//myfile.open ("..\\Output\\cells_11_13_1.txt");
	myfile.open(fn);
	myfile << istream.str();
	myfile.close();
}

/** 
* Run a simulation based on set of boolean parameters
*/
int RunSimulation(Combo c, int it){

	//Set up hallmarks for analysis
	//They are stored in the combo I passed in
	//The code for the combination of hallmarks is as well
	//The code is used in the write to file function
	HALL1 = c.h1;
	HALL2 = c.h2;
	HALL3 = c.h3;
	HALL4 = c.h4;
	HALL5 = c.h5;
	//HALL6 = c.h6;
	//HALL7 = c.h7;
	//HALL8 = c.h8;
	ENABLE2 = c.e2;
	EMERGE1 = c.em1;
	EMERGE2 = c.em2;

	//Create the grid the cells are on
	std::vector<std::vector<Cell *> > grid;
	grid.resize(GRID_SIDE);
	for (int i = 0; i< GRID_SIDE; i++){
		grid[i].resize(GRID_SIDE);
	}
	for (int i=0; i<GRID_SIDE; i++) {
		for (int j=0; j<GRID_SIDE; j++){
			grid[i][j] = NULL;
		}
	}

	////Set up the oxygen grid
	BinaryFluid binaryGrid = BinaryFluid(GRID_SIDE);
	//binaryGrid.updateBinaryFluid();


	//Use a vector to keep track of all of the cells.. just add them onto/into the vector when created
	list<Cell*> * allCells = new list<Cell*>();

	//Keep track of total number of cells and cancer cell
	//End simulation if END_PERCENT of cells are cancerous
	int cellCount = 0;
	int cancerCount = 0;
	int aliveNormal = 1;
	double cancerPercent = 0;

	//Set time to 0
	int time = 0;
	//Only update binary fluid if the time has changed
	bool timeChanged = false;

	//Priority queue 
	priority_queue<Event, vector<Event>, CompareEvent> events;
	//Create initial cell and add pointer to grid
	Cell *start = new Cell();
	//Add cell to the grid
	grid[start->geti()][start->getj()] = start;

	//Add cell to the list and increment the total cell count
	allCells->push_back(start); cellCount++;
	Event *firstEvent = new Event(*start, time);
	events.push(*firstEvent); //Add first event to queue

	Event currentEvent; //is it okay these are on the stack? I'm reassigning them each time through loop..
	Cell * currentCell; //1/24
	int currentTime;
	Cell * daughterCell;
	Event * newEvent1;
	Event * newEvent2;

	//temp int counter
	int counter = 0;
	
	//Main body of simulator
	//While there are still events to process, the simulation hasn't been going for too long and the tumour isn't completely cancerous, keep simulating
	//while ((!(events.empty())) && (counter<20000) && (cancerPercent < END_PERCENT)){
	while ((!(events.empty())) && (counter<40000)){
		//Get next event
		currentEvent = events.top();
		//Get the current cell
		currentCell = currentEvent.getCell();
		//Only do anythign with it if the cell is alive (dead cells don't grow or consume oxygen)
		//if (currentCell->getState() == ALIVE){
		if (currentCell->isAlive()) {
			//Get the time of this operation from the queue
			timeChanged = false;
			currentTime = currentEvent.getTime();
			if (currentTime > time){
				time = currentTime;
				timeChanged = true;
			}

			//Check if it died (random cell death)
			if (currentCell->died()) { //if it dies, move on to next cell.. 
				events.pop();
				counter++;
				continue;
			}
			//Check for apoptosis
			if (currentCell->isMutated()) {
				if (currentCell->apoptosis()==true) { //if it dies via apoptosis, move on
					events.pop(); //This removes the older event
					counter++;
					continue;
				}
			}

			//Next, do the 3 mitosis checks:
			//1) Check if within growth factor check OR if hallmark is on
			//Hallmark 1 advantage
			bool canGrow = false;
			if (currentCell->selfGrows() || currentCell->withinGrowthFactorRange()){
				canGrow = true;
			}
			//2) Check if there is space or if it ignores it 
			bool space = false;
			//The hasSpace() and compete() functions return -1 if no space
			//or return a direction for movement if there is space
			//Hallmark 3 advantage is in here
			Direction newDir = (Direction)currentCell->hasSpace();
			if (newDir != -1){
				space = true;
			}
			//Check if the telomere is still around
			//Hallmark 4 advantage
			bool telomere = false;
			if (currentCell->getTelomere() > 0 || (currentCell->ignoresTelomere()==true)){
				telomere = true;
			}
			else {
				//If the cell is out of telomere then it is dead
				currentCell->setToDead();
				events.pop(); //This removes the older event
				counter++;
				continue;
			}
			//Check Hallmark 5 .. within blood range or has angiogenesis turned on
			bool blood = false;
			
			if (currentCell->withinBloodRange() == true){
				blood = true;
			}
			//else {
			//	//If cell is out of blood (and therefore oxygen) it dies
			//	currentCell->setToDead();
			//	events.pop(); //This removes the older event
			//	counter++;
			//	continue;
			//}
			//Apply immune system
			if (currentCell->killedByImmune()){
				events.pop(); //This removes the older event
				counter++;
				continue;
			}



			//2/3
			//Check if the cell has oxygen
			//Check that the cell is still alive Jenna to-do
			double oxygenAmount = binaryGrid.getOxygenValue(currentCell->geti(), currentCell->getj());
			bool enoughOxygen = currentCell->checkOxygen(oxygenAmount);
			//if there is enough oxygen, consume it
			if (enoughOxygen){
				binaryGrid.consumeOxygen(currentCell->geti(), currentCell->getj(), currentCell->getRequiredOxygen());
				currentCell->markConsumedOxy(true);
				//I think marking it this way but *not* checking it here 
				//eans it will consume each time it is actively involved in
				//a mitosis and comes up in the queue, but won't re-consume
				//when we are doing the big oxygen check.. Mark?
				//Perhaps this should be inside? Or just last once we know
				//it is making it through all the changes?
			}
			//If there wasn't enough oxygen the cell should be quiescent (if enough for that)
			//Then the cell shoudl break out of this loop
			else if (oxygenAmount >= OXYGEN_CON_QUIS) {
				currentCell->setState(QUIS);
				binaryGrid.consumeOxygen(currentCell->geti(), currentCell->getj(), OXYGEN_CON_QUIS);
				currentCell->markConsumedOxy(true);
				//Don't go forward with mitosis if not enough oxygen
				events.pop(); //This removes the older event
				counter++;
				continue;
			}
			//Not even enough oxygen for quiescense..
			else {
				//JENNA TO-DO: check that this changes the cell state and the alive bool
				currentCell->setState(NEC); //cell dies via apoptosis
				events.pop(); //This removes the older event
				counter++;
				continue;
			}

			//Check if the cell is "trapped".. if it has no space,
			//is outside of growth factor and outside of angiogenesis
			//and _not_ cancerous, it is trapped and should die
			if (currentCell->checkTrapped()){
				events.pop();
				counter++;
				continue;
			}
			//Lastly, make sure cell is alive
			bool cellAlive = currentCell->isAlive();
			//If all of the above conditions are satisfied, we perform mitosis 
			if (canGrow && space && telomere && blood && cellAlive && enoughOxygen){

				//Perform mitosis
				//Create a daughter cell and add it to the list; increment counter
				daughterCell = new Cell(*currentCell); //Call copy constructor BUT 
				//set its own i,j, telomere and neighbours later on
				//TO-DO Jenna updated telomere?

				//This daugther cell needs to be a copy... except for the telomere.. it should
				//start out high again...
				//Make sure the daughter cell is getting the parent cell parameters like
				//any cancer mutations

				//set the i,j,k of the new cell
				daughterCell->calcIJK(currentCell->geti(), currentCell->getj(), currentCell->getk(), newDir);
				//Add new cell to grid
				grid[daughterCell->geti()][daughterCell->getj()] = daughterCell;
				//Update cells neighbour list
				daughterCell->updateOwnNeighbours(grid);

				//Update neighbours
				//Current cell needs to be udpated with the daughter cell information.. that will happen in update all
				currentCell->updatedNeighbour(newDir, daughterCell); 
				daughterCell->updateOtherNeighbours(grid);

				//Chance for mutation and decrease telomere
				currentCell->decreaseTelomere();
				//daughterCell->setTelomere(100); //New cells start with full length telomere
				daughterCell->mitosisOccured();
				currentCell->mitosisOccured();


				//push this cell onto the list of all cells, and increase cell count
				allCells->push_back(daughterCell); cellCount++;
				//Create the new event and then push it onto the queue
				newEvent1 = new Event(*daughterCell, time, newDir); //is this okay? Mark? Because I reassign this pointer and put it on the queue... 
				//will I be losing references to that which is on the queue?
				events.push(*newEvent1);
				//Jenna added this if
				if (currentCell->hasAnySpace() || currentCell->ignoresGrowthInhibition()==true) {
					newEvent2 = new Event(*currentCell, time);
					events.push(*newEvent2);
				}

			}
		}

		events.pop(); //This removes the older event
		counter++;

		//Check with Mark as to how many steps
		//Also counts how many cancer cells there are
		if ((timeChanged == true) && (time % 25 == 0)){
			
			std::list<Cell *>::const_iterator iterator = allCells->begin();
			std::list<Cell *>::const_iterator iterator2 = allCells->end();
			double oxygenAmount;
			bool enoughOxygen;
			//Go through and have all cells consume oxygen
			//If the cell is alive, it should consume healthy cell oxygen levels
			//If it is quiescent it should consume half that
			//If its trapped, it dies
			for (; iterator != iterator2; ++iterator) {
				//TO-DO check that my isAlive check here is working correctly
				//All living cells should consume some level of oxygen
				if (((*iterator)->isAlive()) && ((*iterator)->hasAlreadyConsumedOxy() == false)) {
					oxygenAmount = binaryGrid.getOxygenValue((*iterator)->geti(), (*iterator)->getj());
					enoughOxygen = (*iterator)->checkOxygen(oxygenAmount);
					//If the cell has enough oxygen, conusme it
					if (enoughOxygen){
						binaryGrid.consumeOxygen((*iterator)->geti(), (*iterator)->getj(), (*iterator)->getRequiredOxygen());
						(*iterator)->markConsumedOxy(true);
					}
					//Not enough oxygen to consume, so the cell becomes quiescent
					else if (oxygenAmount >= OXYGEN_CON_QUIS) {
						binaryGrid.consumeOxygen((*iterator)->geti(), (*iterator)->getj(), OXYGEN_CON_QUIS);
						(*iterator)->setState(QUIS);
					}
					//Not enough oxygen to survive.. cell becomes necrotic
					else {
						(*iterator)->setState(NEC);
					}
				} //set to false hre so we don't have to go through again //check it with mark
				//Reset it for next round
				(*iterator)->markConsumedOxy(false);
				//Count number of alive normal cells, and alive cancer cells
				//Alive is 0
				if ((*iterator)->isMutated() && ((*iterator)->getState() == 0)){
					cancerCount++;
				}
				if (((*iterator)->isMutated() == false) && ((*iterator)->getState() == 0)){
					aliveNormal++;
				}
			}
			cancerPercent = cancerCount/aliveNormal; //TO-DO change this to be ALIVE cancer and ALIVE regular... hm...
			//if (cancerPercent > END_PERCENT){
				//return 1;
			//}
			cancerCount = 0; //Reset it because we don't want to be adding to it
			aliveNormal = 1; //keep at 1 so we don't get divide by 0 error if there are no alive cells left

			//Update binary fluid
			//2/10

			binaryGrid.updateBinaryFluidSmall();
		}

	//	Write out every x timesteps
		if (counter % 500 == 0){
			std::list<Cell *>::const_iterator iterator = allCells->begin();
			std::list<Cell *>::const_iterator iterator2 = allCells->end();
			WriteToFile(iterator, counter, iterator2, c.code, it, time);
			cout << time <<"\n";
		}

	}

	std::list<Cell *>::const_iterator iterator = allCells->begin();
	std::list<Cell *>::const_iterator iterator2 = allCells->end();
	WriteToFile(iterator, counter, iterator2, c.code, it, time);

	for (std::list<Cell *>::const_iterator iterator = allCells->begin(), end = allCells->end(); iterator != end; ++iterator) {
		//(*iterator)->print();
		delete *iterator;
	}

	delete allCells;
	//delete start;
	//delete firstEvent;

	return 0;
}

list<Combo> * allCombos = new list<Combo>();

void generateGroups() {
	Combo c1 = {"a", 1, 1, 1, 1, 1, 1, 1, 1, 1};	allCombos->push_back(c1);
	//Combo c2 = {"b", 0, 0, 1, 1, 1, 1, 1, 1};	allCombos->push_back(c2);
	//Combo c3 = {"c", 0, 1, 0, 1, 1};	allCombos->push_back(c3);
	//Combo c4 = {"d", 0, 1, 1, 0, 1};	allCombos->push_back(c4);
	//Combo c5 = {"e", 0, 1, 1, 1, 0};	allCombos->push_back(c5);


	//Combo c6 = {"f", 1, 0, 0, 1, 1};	allCombos->push_back(c6);
	//Combo c7 = {"g", 1, 0, 1, 0, 1};	allCombos->push_back(c7);
	//Combo c8 = {"h", 1, 0, 1, 1, 0};	allCombos->push_back(c8);

	//Combo c9 = {"i", 1, 1, 0, 0, 1};	allCombos->push_back(c9);
	//Combo c10 = {"j", 1, 1, 0, 1, 0};	allCombos->push_back(c10);

	
	//Combo c11 = {"k", 1, 1, 1, 0, 0};	allCombos->push_back(c11);

}



/** 
* Main
* Creates the groups of hallmarks
* Loops through each and runs simulation for each set of hallmarks
*/
int _tmain(int argc, _TCHAR* argv[])
{

	srand(RAND_SEED);

	

	//Turn on loop for multiple comparissions 
	//for (int i=0; i<10; i++){
	generateGroups();
	int howEnd = 0;
	while (allCombos->empty() != true){
		Combo c = allCombos->front();
		howEnd = RunSimulation(c, 0);
		allCombos->pop_front();
	}
	
	//Now generate the pictures and such
	std::string test = "python analyzeDataApril20.py";
	system(test.c_str());

	WriteParams(howEnd);
}


//Cell *x = new Cell();
//x->print();

//Cell x2;
//x2.setMutationRate(10);

//CancerCell * y = new CancerCell();
//y->print();
//printf("\n %d", y->getMutationRate());


//x->setMutationRate(5);
//x->setTelomere(10);

//Event test(*x, 10);
//test.setTime(5);
//Cell t = test.getCell();
//printf("cells mutation rate: %d", t.getMutationRate());

////Here I am attempting to make a Cell pointer and point it to a cancer cell
////I have avoidsApoptosis as virtual in my Cell class
////Then it is implemented in teh cancer cell
////I thought I could then do xy->avoidsApop but it is giving me a memory error
////This doesn't work.. when I do cancer1.avoidsApoptosis it's like cancer1 doesn't exist...
////CancerCell cancer1;
////Cell * xy = &cancer1;
////xy->setTelomere(5);
////cancer1.setAvoidApop(true); //works
////xy->avoidsApoptosis(); //memory error
////cancer1.avoidsApoptosis();
////printf("avoids apoptosis: %s", xy->avoidsApoptosis());
////Ask Mark about the above issue...


//int temp[] = {1,1,1,1,1,1,1,1};
//y->setNeighbours(temp, 8);

//const int * neigh = y->getNeighbours(); //Here testing stuff.... I think the get neighbours header and implementation is wrong
//printf("%d", neigh[0]);

////This creates a 2D vector of health cells... Not sure if we want to do that though... 
////vector<vector<Cell> > cells(100);
////for (int i=0; i<100; i++){
////	cells[i].resize(100);
////}

///*


//Mark
//CancerCell cancer1;
//Cell * xy = &cancer1;
//xy->setTelomere(5);
//cancer1.setAvoidApop(true); //works
//bool testbo = xy->avoidsApoptosis(); //why does "this" look funny? Mark
//cancer1.avoidsApoptosis();
//printf("avoids apoptosis: %d", testbo);

//Testing the neighbours
//Mark does this seem right?
//Cell test;
//test.setTelomere(5);
//Cell * neigh = test.getNeighbours();
//for (int i=0; i<8; i++){
//	if (neigh)
//		printf("huh");
//	else
//		neigh++;
//}

////Test space
//Cell test;
//Cell north;
//Cell south;
//Cell southWest;
//Cell * neigh[8] = {&north, NULL, &south, NULL, NULL, NULL, &southWest, NULL};
//test.setNeighbours(neigh);
//int x = test.hasSpace();

//Priority Queue Testing


//Test neighbour idea
//Cell cell1;
//Cell cell2;
//cell2.calcIJK(cell1.geti(), cell1.getj(), cell1.getk(), 2);
//Cell cell3;
//cell3.calcIJK(cell2.geti(), cell2.getj(), cell2.getk(), 1);
//Cell * temp = &cell3;
//cell1.updatedNeighbour((Direction)5, temp);


//Attempting to do a grid with pointers to pointers
////Get a grid of pointers for the cells
//	Cell*** cellGrid = new Cell**[GRID_SIDE];
//	for(int i = 0; i < GRID_SIDE; i++)
//		cellGrid[i] = new Cell*[GRID_SIDE];
//
//	for (int i=0; i<GRID_SIDE; i++){
//		for (int j=0; j<GRID_SIDE; j++){
//			cellGrid[i][j] = NULL;
//		}
//	}
//	//Get a pointer to the grid
//	Cell (*p_cellGrid) = **cellGrid;
//in cell.h:
//int updateOwnNeighbours(Cell (*cellGrid)[GRID_SIDE][GRID_SIDE]);

//in cell
/**
Use the pointer grid to update own neighbours list, then you can call updateAllNeighbours
*/
//int Cell::updateOwnNeighbours(Cell (*cellGrid)[GRID_SIDE][GRID_SIDE]){
//	for (int i=0; i<NUM_NEIGHBOURS; i++){
//		if (cellGrid[i][j-1] != NULL)
//			neighbours[0] = cellGrid[i][j-1];
//		if (cellGrid[i+1][j] != NULL)
//			neighbours[1] = cellGrid[i+1][j];
//		if (cellGrid[i][j+1] != NULL)
//			neighbours[2] = cellGrid[i][j+1];
//		if (cellGrid[i-1][j] != NULL)
//			neighbours[3] = cellGrid[i-1][j];
//		if (cellGrid[i+1][j-1] != NULL)
//			neighbours[4] = cellGrid[i+1][j-1];
//		if (cellGrid[i+1][j+1] != NULL)
//			neighbours[5] = cellGrid[i+1][j+1];
//		if (cellGrid[i-1][j+1] != NULL)
//			neighbours[6] = cellGrid[i-1][j+1];
//		if (cellGrid[i-1][j-1] != NULL)
//			neighbours[7] = cellGrid[i-1][j-1];
//	}
//	return 0; //TO-DO error checking
//}

////Testing cancer cell properties for normal cell
//Cell *x = new Cell();
//x->print();

//Cell x2;
//x2.setMutationRate(10);

//x->setIgnoreGrowth(true);
//printf("Ignores: %d", x->ignoresGrowthInhibition());
//x->setSelfGrowth(true);
//printf("self growth? %d", x->selfGrows());

//I put this into the space code...
//else if (currentCell.ignoresGrowthInhibition()){
//	newDir = (Direction)currentCell.compete();
//	if (newDir != -1)
//		space = true;
//}
//3) Check telomere

//Write to file code:
////Get all i,j coordinates
//stringstream istream;
////11/04
////self grows, ignores inhibition, avoids apoptosis, ignores telo, deregulates, avoid immunity, genomic instability
//int muts[] = {0, 0, 0, 0, 0, 0, 0};
//int cellTotal = 0;
//istream << "i, j, state, mut, selfGrow, ignoresGrow, avoidsAp, ignoresT, dereg, avoidsI, genomeUn \n";
////Iterate over the cells and get list to draw
//for (std::list<Cell *>::const_iterator iterator = allCells->begin(), end = allCells->end(); iterator != end; ++iterator) {
//	//(*iterator)->print();
//	istream << (*iterator)->geti() << " " << (*iterator)->getj() << " " << (*iterator)->getState() << " " << (*iterator)->isMutated() << " ";
//	istream << (*iterator)->selfGrows() << " " << (*iterator)->ignoresGrowthInhibition() << " " << (*iterator)->avoidsApoptosis() << " ";
//	istream << (*iterator)->ignoresTelomere() << " " << (*iterator)->deregulatesCellEnergetics() << " " << (*iterator)->avoidsImmunity() << " " << (*iterator)->genomeUnstable() <<"\n";
//}

////Write the cells to a file
//istream << "\n";
//ofstream myfile;
//myfile.open ("..\\Output\\cells_11_11_1.txt");
//myfile << istream.str();
//myfile.close();

//return 0;