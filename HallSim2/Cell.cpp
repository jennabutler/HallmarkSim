/*
* Cell.cpp
*
*  Created on: Oct 1, 2013
*      Author: jenna
*/


#include "stdafx.h"
#include "Cell.h"
#include "Global.h"
#include "Constants.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

namespace std {

	enum Direction {
		N, E, S, W, NE, SE, SW, NW
	};

	enum State {
		ALIVE, DEAD, QUIS, APOP
	};

	/** 
	\class Cell
	\brief The cell class
	*/

	Cell::Cell() {
		//Telomere, mutation and death rates from Santos 2012
		telomere = 100;
		evadeApoptosis = 10;
		mutationRate = MUTATION_RATE;
		geneticInstabilityFactor = 5;
		randomDeathRate = 10000;
		competition = 10; 
		//Inialize neighbour pointers to be null
		for (int c = 0; c < NUM_NEIGHBOURS; c++){
			neighbours[c] = NULL;
		}

		mutated = false;
		state = ALIVE;

		//Cell start at the center of the grid
		i = (GRID_SIDE/2)+1;
		j = (GRID_SIDE/2)+1;
		if (DIM == 3) 
			k = (GRID_SIDE/2)+1;
		else
			k = -1;

		//Cancer Cell Properties... only turned on once it is cancerous
		selfGrowth = false;
		ignoreGrowthInhibition = false;
		apopLiklihood = 0;
		avoidApoptosis = false;
		ignoreTelomere = false;
		susAngio = false;
		//deregulateCellEnergetics = false;
		immuneDeathLiklihood = 0;
		avoidImmunity = false;
		genomicInstability = false;

		//Each round there is a flag for whether or not a cell has
		//consumed any oxygen...
		//Each time a cell mutates it consumes oxygen
		//Then once every few time steps it consumes
		//If it has already consumed because of mutation
		//it doesn't do it in the time step check
		consumed = false;

	}

	Cell::~Cell()
	{
		// TODO Auto-generated destructor stub
	}

	/**
	* @return int mutation rate
	*/
	int Cell::getMutationRate() const
	{
		return mutationRate;
	}

	/**
	* Calculates whether a cell dies
	* Cells have chance 1/randomDeathRate at dying each time through their lifecycle
	* @return bool whether the cell died or not
	*/
	bool Cell::died() {
		int x = rand() % randomDeathRate + 1;
		if (x == randomDeathRate) {
			//state = DEAD;
			setToDead();
			return true;
		}
		return false;
		//If a cell dies, it stays where it is.
		//I may change it so that if it dies that space becomes available
		//If I do that, I need to add some code to handle updatin the neighbours
	}

	/**
	* Apoptosis
	* A cell with n hallmarks mutated has an extra n/e liklihood of dying each
	* cell cycle unless the "evade apoptosis" hallmark is on
	*/
	bool Cell::apoptosis(){
		//Hallmark 3 advantage.. never will die via apoptosis
		if (avoidApoptosis == true){
			return false;
		}
		int mutCount = 0;
		if (selfGrowth) mutCount++;
		if (ignoreGrowthInhibition) mutCount++;
		if (avoidApoptosis) mutCount++;
		if (ignoreTelomere) mutCount++;
		//if (deregulateCellEnergetics) mutCount++;
		if (avoidImmunity) mutCount++;
		if (genomicInstability) mutCount++;
		if (susAngio) mutCount++;
		//float x = rand();
		//float prob = (float)(mutCount)/(float)(evadeApoptosis);
		//if (x < prob){
		////if (rand() < (float)(mutCount)/(float)(evadeApoptosis)) {
		//	//apoptosis
		//	state = DEAD;
		//	return true;
		//}
		double apopChance = ((float)(mutCount)/evadeApoptosis)*evadeApoptosis;
		int x = rand() % evadeApoptosis + 1;
		if (x <= apopChance){ //Cell dies via apoptosis
			//state = DEAD;
			setToDead();
			return true;
		}
		return false;
	}


	/**
	* Accessors
	*/

	/**
	* @return the randomDeathRate (int)
	*/
	int Cell::getRandomDeathRate() const
	{
		return randomDeathRate;
	}

	/**
	* @return the current telomere length (int)
	*/
	int Cell::getTelomere() const
	{
		return telomere;
	}

	/**
	* Helper method to get the neighbour list
	* @return a pointer to the first neighbour
	*/
	Cell * Cell::getNeighbours()
	{
		Cell * neigh = neighbours[0];
		return neigh; //returns a pointer to the first item in the neighbour list
	}

	/**
	* @return current state of the cell (an int, but an enum)
	*/
	int Cell::getState() const
	{
		return state;
	}

	/**
	* @return whether the cell is mutated or not.. normal cells are not, cancer cells are (bool)
	*/
	bool Cell::isMutated() const
	{
		return mutated;
	}

	/**
	* Sets the mutation rate for the cell
	* @param new rate of mutation for the cell
	*/
	void Cell::setMutationRate(int mutationRate)
	{
		this->mutationRate = mutationRate;
	}

	/**
	* Sets the random death rate
	* @param int for the new random death rate
	*/
	void Cell::setRandomDeathRate(int randomDeathRate)
	{
		this->randomDeathRate = randomDeathRate;
	}

	/**
	* Sets the telomere length (starts at 100)
	* @param new telomere length
	*/ 
	void Cell::setTelomere(int telomere)
	{
		this->telomere = telomere;
	}

	/**
	* Sets the neighbour list
	* @param an array of pointers to neighbouring Cell objects
	*/
	void Cell::setNeighbours(Cell * newNeighbours[])
	{
		for (int i=0; i<NUM_NEIGHBOURS; i++)
		{
			neighbours[i] = newNeighbours[i];
		}
	}

	/**
	* Prints out details on the cell
	*/
	void Cell::print()
	{
		printf("Cell at i: %d, j: %d, mutation:  %d \n", i, j, mutated);
		//printf("%s", "in parent cell");
	}



	/**
	* Method to calculate if the cell is within the range of the
	* growth factor
	* TO-DO Change it to a circle or something, not a square
	* @return bool true if within range, false otherwise
	* TO-DO: They used 95% of inner space but that doesn't make sense to me
	*/ 
	bool Cell::withinGrowthFactorRange(){
		//Hallmark 1 preference
		if (selfGrowth == true)
			return true;
		//Equation of a circle is (x-h)^2 + (y-k)^2 = r2
		//for a circle centered at h,k with radius r
		//This makes a circle of growth factor
		int within = rand() % 10 + 1; //drop off the chances of it surviving
		if (((pow((i-GROWTH_H),2)) + (pow((j-GROWTH_K),2))) <= (pow(GROWTH_RADIUS,2))){
			return true;
		}
		else if (((pow((i-GROWTH_H),2)) + (pow((j-GROWTH_K),2))) <= (pow(GROWTH_RAIDUS_2,2))) {
			if (within < 6)
				return true;
			return false;
		}
		else if (((pow((i-GROWTH_H),2)) + (pow((j-GROWTH_K),2))) <= (pow(GROWTH_RADIUS_3,2))) {
			if (within < 4)
				return true;
			return false;
		}
		else if (((pow((i-GROWTH_H),2)) + (pow((j-GROWTH_K),2))) <= (pow(GROWTH_RADIUS_4,2))) {
			if (within < 2)
				return true;
			return false;
		}
		else {
			return false;
		}
		//Make parabolic?
	}

	/**
	* Method to calculate if the cell is within the range of the
	* blood vessels for oxygen
	* @return bool true if within range, false otherwise
	* TO-DO: If not within blood range, should it die?  handled in sim
	*/ 
	//Should be a parabola?
	bool Cell::withinBloodRange(){
		//Hallmark 5 check
		//If the cell has sustained angiogenesis on, it doesn't matter if it is within the range of blood
		if (susAngio == true) {
			return true;
		}
		//If the cell's neighbours have sustained angiogenesis, it gets blood regardless
		else if (neighbourSusAngio() == true) {
			return true;
		}
		int within = rand() % 10 + 1; //drop off the chances of it surviving
		if (((pow((i-BLOOD_H),2)) + (pow((j-BLOOD_K),2))) <= (pow(BLOOD_RADIUS,2))){
			return true;
		}
		else if (((pow((i-BLOOD_H),2)) + (pow((j-BLOOD_K),2))) <= (pow(BLOOD_RADIUS_2,2))) {
			if (within < 6)
				return true;
			return false;
		}
		else if (((pow((i-BLOOD_H),2)) + (pow((j-BLOOD_K),2))) <= (pow(BLOOD_RADIUS_3,2))) {
			if (within < 4)
				return true;
			return false;
		}
		else if (((pow((i-BLOOD_H),2)) + (pow((j-BLOOD_K),2))) <= (pow(BLOOD_RADIUS_4,2))) {
			if (within < 2)
				return true;
			return false;
		}
		else {
			return false;
		}
	}


	
	/**
	* Method to determine if a cancer cell is killed off by the immune system
	* Cell must be mutated to be killed
	* If cell is on an angiogenesis line, it has a higher chance of being killed off
	* All cancer cells have some chance, but those on a blood line have a higher chance
	* A cells individual chance is resprested by 
	*/bool Cell::killedByImmune(){
		if (mutated == false) {
			//Not a cancer cell, so will not be killed by immune system
			return false;
		}
		//Next check if it has anigogenesis turned on, or its neighbours do
		else {
			//So the cell is on a blood path and is cancerous
			//Generate a number between 1 and the likelihood
			int immuneDeathChance = rand() % immuneDeathLiklihood + 1;
			if (immuneDeathChance == immuneDeathLiklihood){
				setToDead(); //TO-DO
				return true; //Cell has been killed
			}
			else {
				return false;
			}
		}
	}

	/**
	* Method to calculate if a cell is close to another cell
	* who is alive and with angiogenesis activated
	* @return bool true if within range, false otherwise
	*/ 
	bool Cell::neighbourSusAngio(){
		//Need to check if any neighbours around me have that... Hmm...
		if (neighbours[0] != NULL) {
			if (neighbours[0]->angiogenesis()) {
				return true;
			}
			////Now get you neighbours neighbours and check if they are angiogenic
			////Mark... still an error. It does increase by 4, but now that isn't null
			//Cell ** neighbour = (Cell **)(neighbours[0])->getNeighbours();
			//
			//for (int i=0; i<8; i++){
			//	if (neighbour != NULL) {
			//		if ((*neighbour)->angiogenesis())
			//			return true;
			//	}
			//	neighbour = neighbour++; 
			//}
		}
		if (neighbours[1] != NULL) {
			if (neighbours[1]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[2] != NULL) {
			if (neighbours[2]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[3] != NULL) {
			if (neighbours[3]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[4] != NULL) {
			if (neighbours[4]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[5] != NULL) {
			if (neighbours[5]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[6] != NULL) {
			if (neighbours[6]->angiogenesis()) {
				return true;
			} 
		}
		if (neighbours[7] != NULL) {
			if (neighbours[7]->angiogenesis()) {
				return true;
			} 
		}
		return false;
	}

	/**
	* Determines if the cell has space to grow (a null pointer in neighbour list)
	* Also decides where the cell should grow
	* TO-DO: Pointers on the boundaries need to be null.. somehow need to implement that
	* @return the int location where the cell should grow into
	*/ 
	int Cell::hasSpace(){
		//srand(1234); //11/04
		//bool spaces[NUM_NEIGHBOURS];
		bool space = false;
		for (int i=0; i<NUM_NEIGHBOURS; i++){
			if ((neighbours[i])==NULL) {
				space = true;
			}
		}
		//If there are no spots,check for competition, or return -1
		if (space == false) {
			//Hallmark 2 Preference
			if (ignoreGrowthInhibition == true) {
				return compete(); //Compete for a space if no space and you have this mutation
			}
			return -1;
		}

		//If there is space,
		//Go through a pick a random one of the true spots...
		bool foundSpot = false;
		int newSpot=0;
 		while (!(foundSpot)){
			//printNeighbours();
			newSpot = rand() % NUM_NEIGHBOURS + 0;
			//Make sure the newSpot is within the range of our neighbour list
			//so we don't accidentally index into random memory
			if ((newSpot > NUM_NEIGHBOURS-1) || (newSpot < 0))
				return -2; //problem with code;
			if (neighbours[newSpot] == NULL) {
				return newSpot; //newspot is an index in the list... it should also be a enum for direction
			}
		}
		//If you can't find a spot for some reason, return -1
		return -1;
	}

	//Jenna
	bool Cell::hasAnySpace(){
		bool space = false;
		for (int i=0; i<NUM_NEIGHBOURS; i++){
			if (neighbours[i] == NULL)
				space = true;
		}
		return space;
	}

	/**
	* When a new cell is created, there is a chance for mutation
	* Each hallmark has a 1/m chance of mutation, where m is the base mutation rate
	*JENNA TO-DO chance the comparrision to be equal to mutation rate and make sur
	*/
	void Cell::mitosisOccured(){
		int mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate) && (HALL1 == true)){
			//Turn on hallmark 1
			selfGrowth = true;
			mutated = true;
			updateImmuneDeathLiklihood();
		}
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate) && (HALL2 == true)){
			//Turn on hallmark 2
			ignoreGrowthInhibition = true;
			mutated = true;
			updateImmuneDeathLiklihood();

		}
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate) && (HALL3 == true)){
			//Turn on hallmark 3
			avoidApoptosis = true;
			mutated = true;
			updateImmuneDeathLiklihood();
		}
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate)&& (HALL4 == true)){
			//Turn on hallmark 4
			ignoreTelomere = true;
			mutated = true;
			updateImmuneDeathLiklihood();
		}
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate) && (HALL5 == true)){
			//Turn on hallmark 5
			susAngio = true;
			mutated = true;
			updateImmuneDeathLiklihood();
		}
		//TO-DO Look up a proper parameter for this value (Abbott?)
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == (mutationRate)) && (EMERGE2 == true)){
			//Turn on emerging hallmark 2
			avoidImmunity = true;
			//Now that mutated, update the liklihood cell will be killed by immune system
			updateImmuneDeathLiklihood();
			mutated = true;
		}
		//No effect in the code 
		//mutHall = rand() % mutahtionRate + 0;
		//if ((mutHall  > (mutationRate-5))&& (HALL8 == true)){
		//	//Turn on hallmark 8
		//	avoidImmunity = true;
		//	mutated = true;
		//}
		mutHall = rand() % mutationRate + 1;
		if ((mutHall == mutationRate)&& (ENABLE2 == true)){
			//Turn on enabling characteristc 2
			//First time it gets turned on, 
			//update hte mutation rate to be faster
			if (genomicInstability == false){
				mutationRateUpdate();
			}
			genomicInstability = true;
			mutated = true;
			updateImmuneDeathLiklihood();
		}
		//Everytime through the loop there is an increased chance of
		//changing the mutation rate if the cell is genetically unstable
		//mutHall = rand() % mutationRate + 0;
		//if (genomicInstability && (mutHall > (mutationRate-5))){
		//	mutationRateUpdate();
		//}
	}

	/**
	* Method to update the immune death liklihood parameter
	* Called whenever a cell mutates
	*/
	void Cell::updateImmuneDeathLiklihood(){
		//If you are on the blood system, you have a higher chance
		//of being killed by the immune system (1 in 100)
		if (susAngio == true || neighbourSusAngio() == true){
			immuneDeathLiklihood = 100;
		}
		//If not on the blood system, a 1 in 1000 chance of death by immune system
		else {
			immuneDeathLiklihood = 1000;
		}
	}


	/**
	* Method to check if a cell is trapped
	* If a normal cell is at the center and has no space,
	* no blood and no growth factor it dies
	* TO-DO: Discuss with Mark if this makes sense
	*/
	void Cell::checkTrapped(){
		//Should we do the blood one? Mark
		//if (!(mutated) && !(withinBloodRange()) && (!(withinGrowthFactorRange())) && (hasSpace() == -1)){
		if ((!(mutated)) && (!(withinGrowthFactorRange())) && (hasSpace() == -1)){
			setToDead();
		}
	}

	/** Method to test if enough oxygen to live
	//2/3
	*/
	bool Cell::checkOxygen(double oxygenAmount){
		if (oxygenAmount >= OXYGEN_CONSUMPTION){
			return true;
		}
		//Mark: I'm not sure about this
		//I only want to consume oxygen if there is enough there to remove 
		//from the lbm system
		//If the cell is angiogenic, should it really be removing oxygen 
		//from the system? It should in reality, but our system 
		//doesn't infuse oxygen into the system, so I feel like the 
		//angiogenic ones shouldn't remove it perhaps?
		else if (angiogenesis() == true){
			return true;
		}
		else if (neighbourSusAngio() == true){
			return true;
		}
		else {
			return false;
		}
	}

	/**
	* If genomic instability is turned on, there is a chance to increase the ovearll
	* mutation rate at each division
	*/
	//Enabling characteristic advantage
	void Cell::mutationRateUpdate(){
		//TO-DO how much do we lower m by?
		//int changeAmount = mutationRate/100;
		//mutationRate = mutationRate - changeAmount; //11/04
		mutationRate = mutationRate/geneticInstabilityFactor;
	}

	/**
	* Set the state to dead if telomere runs out
	*/
	void Cell::setToDead(){
		state = DEAD;
	}



	/** 
	* Calculate the new i coordinate for a cell
	* You have the other cells i,j and k coordinate
	* You also have the direction YOU are in.. so if it is north, you are north of the cell you have the coordinates of
	*/
	int Cell::calcIJK(int oldI, int oldJ, int oldK, int direction){
		int newI;
		int newJ;
		int newK = 0; //TO-DO: Implement for 3D
		if (direction==N){
			newI = oldI;
			newJ = oldJ -1;
		}
		else if (direction == E){
			newI = oldI + 1;
			newJ = oldJ;
		}
		else if (direction == S){
			newI = oldI;
			newJ = oldJ + 1;
		}
		else if (direction == W){
			newI = oldI - 1;
			newJ = oldJ;
		}
		else if (direction == NE){
			newI = oldI + 1;
			newJ = oldJ - 1;
		}
		else if (direction == SE){
			newI = oldI + 1;
			newJ = oldJ + 1;
		}
		else if (direction == SW){
			newI = oldI - 1;
			newJ = oldJ + 1;
		}
		else if (direction == NW){
			newI = oldI - 1;
			newJ = oldJ - 1;
		}
		else {
			//If not a normal direction must be an error
			i = -1;
			j = -1;
			return -1;
		}
		//If out of the bounds, return -1 for an error
		if (newI < 0 || newI > GRID_SIDE || newJ < -0 || newJ > GRID_SIDE){
			i = -1;
			j = -1;
			return -3; //hitting the edge
		}
		i = newI;
		j = newJ;
		return 0;
	}

	/**
	* Update neighbour function
	* This is called by a new cell that is going through its neighbours
	* If it has a neighbour, it tells that cell that it is here
	* @param neighDir The direction is the direction the NEW CELL is in relatively to the cell it is contacting
	* @return -1 if error (outside bounds)
	*/
	int Cell::updatedNeighbour(int neighDir, Cell * newNeigh){
		neighbours[neighDir] = newNeigh;
		return 0; //make this -1 if error TO-DO
	}

	/**Update its own neighbour list after creation **/
	int Cell::updateOwnNeighbours(vector< vector<Cell *> > &grid){
		//First erase any that were accidentally populated
		//TO-DO JENNA test this
		neighbours[0] = grid[i][j-1];
		neighbours[1] = grid[i+1][j];
		neighbours[2] = grid[i][j+1];
		neighbours[3] = grid[i-1][j];
		neighbours[4] = grid[i+1][j-1];
		neighbours[5] = grid[i+1][j+1];
		neighbours[6] = grid[i-1][j+1];
		neighbours[7] = grid[i-1][j-1];

		/*if (grid[i][j-1] == NULL) 
			neighbours[0] = NULL;
		if (grid[i+1][j] == NULL)
			neighbours[1] =NULL;
		if (grid[i][j+1] == NULL)
			neighbours[2] =NULL;
		if (grid[i-1][j] == NULL)
			neighbours[3] = NULL;
		if (grid[i+1][j-1] == NULL)
			neighbours[4] = NULL;
		if (grid[i+1][j+1] == NULL)
			neighbours[5] = NULL;
		if (grid[i-1][j+1] == NULL)
			neighbours[6] = NULL;
		if (grid[i-1][j-1] == NULL)
			neighbours[7] = NULL;

		if (grid[i][j-1] != NULL) 
			neighbours[0] = grid[i][j-1];
		if (grid[i+1][j] != NULL)
			neighbours[1] = grid[i+1][j];
		if (grid[i][j+1] != NULL)
			neighbours[2] = grid[i][j+1];
		if (grid[i-1][j] != NULL)
			neighbours[3] = grid[i-1][j];
		if (grid[i+1][j-1] != NULL)
			neighbours[4] = grid[i+1][j-1];
		if (grid[i+1][j+1] != NULL)
			neighbours[5] = grid[i+1][j+1];
		if (grid[i-1][j+1] != NULL)
			neighbours[6] = grid[i-1][j+1];
		if (grid[i-1][j-1] != NULL)
			neighbours[7] = grid[i-1][j-1];*/
		return 0;
	}

	/**
	* Go through list of neighbours and call the update function on them
	* TO-DO TEST THIS FUNCTION BIG TIME
	*/
	int Cell::updateOtherNeighbours(vector< vector<Cell *> > &grid){
		//Need a pointer to myself...
		Cell * temp = this;
		if (grid[i][j-1] != NULL) 
			grid[i][j-1]->updatedNeighbour(2, temp);
		if (grid[i+1][j] != NULL)
			grid[i+1][j]->updatedNeighbour(3, temp);
		if (grid[i][j+1] != NULL)
			grid[i][j+1]->updatedNeighbour(0, temp);
		if (grid[i-1][j] != NULL)
			grid[i-1][j]->updatedNeighbour(1, temp);
		if (grid[i+1][j-1] != NULL)
			grid[i+1][j-1]->updatedNeighbour(6, temp);
		if (grid[i+1][j+1] != NULL)
			grid[i+1][j+1]->updatedNeighbour(7, temp);
		if (grid[i-1][j+1] != NULL)
			grid[i-1][j+1]->updatedNeighbour(4, temp);
		if (grid[i-1][j-1] != NULL)
			grid[i-1][j-1]->updatedNeighbour(5, temp);
		return 0; //TO-DO error checking?
	}

	void Cell::printNeighbours(){
		int neigh[8] = {0,0,0,0,0,0,0,0};
		if (neighbours[7] != NULL)
			neigh[7] = 1;
		if (neighbours[0] != NULL)
			neigh[0] = 1;
		if (neighbours[4] != NULL)
			neigh[4] = 1;
		if (neighbours[3] != NULL)
			neigh[3] = 1;
		if (neighbours[1] != NULL)
			neigh[1] = 1;
		if (neighbours[6] != NULL)
			neigh[6] = 1;
		if (neighbours[2] != NULL)
			neigh[2] = 1;
		if (neighbours[5] != NULL)
			neigh[5] = 1;
		//printf("%d, %d, %d \n", neigh[7], neigh[0], neigh[4]);
		//printf("%d, %d, %d \n", neigh[3], 1, neigh[1]);
		//printf("%d, %d, %d \n", neigh[6], neigh[2], neigh[5]);
	}



	/**
	* If no space, a cancer cell will compete with neighbouring cells
	* @return returns the location of the cell it will kill.. -1 if lost the competition
	* TO-DO implement
	* TO-DO make this virtual?
	*/
	int Cell::compete(){
		//has a liklihood of 1/g of killing a neighbour to grow
		int chance = rand() % competition + 1; 
		//Get a random number between 0 and g
		//If the number is g, we will kill a neighbour so return
		//the position of a random neighbours

		if (chance == competition) {
			int newSpot = rand() % NUM_NEIGHBOURS + 0;
			if (newSpot < NUM_NEIGHBOURS && newSpot > -1)
				return newSpot; //new spot is a location in list; acts as a direction
		}
		return -1;
	}

} /* namespace std */