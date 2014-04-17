/*
 * Cell.h
 *
 *  Created on: Oct 1, 2013
 *      Author: jenna
 */
 

#ifndef CELL_H_
#define CELL_H_

#include "Constants.h" 
#include <vector>

namespace std {


class Cell {
public:
    Cell();
	//Copy constructor?
    virtual ~Cell();
    int getMutationRate() const;
    int getRandomDeathRate() const;
    int getTelomere() const;
    int getState() const;
    bool isMutated() const;
    void setMutationRate(int mutationRate);
    void setRandomDeathRate(int randomDeathRate);
    virtual void setTelomere(int telomere);
    Cell * getNeighbours();
    void setNeighbours(Cell * newNeighbours[]);
    bool died();
	bool apoptosis();
	void virtual print();
	int hasSpace();
	int geti() { return i; }
	int getj() { return j; }
	int getk() { return k; }
	int calcIJK(int oldI, int oldJ, int oldK, int direction);
	int updatedNeighbour(int neighDir, Cell * newNeigh);
	//int updateAllNeighbours();
	int updateOwnNeighbours(vector< vector<Cell *> > &grid);
	int updateOtherNeighbours(vector< vector<Cell *> > &grid);
	void printNeighbours();
	bool hasAnySpace(); //jenna
	void mitosisOccured();
	void decreaseTelomere() { telomere--; }
	void mutationRateUpdate();
	void setToDead();
	void setState(int stateType);
	void updateImmuneDeathLiklihood();
	void checkTrapped();

	//Methods not yet implemented correctly
	bool withinGrowthFactorRange();
	bool withinBloodRange();
	bool killedByImmune();
	bool neighbourSusAngio();
	int compete();
	bool checkOxygen(double oxygenAmount, int cellState);
	void markConsumedOxy(bool consumedOxygen) { consumed = consumedOxygen; }
	bool hasAlreadyConsumedOxy() { return consumed; }


	//Virtual methods available in the cancer child
	//virtual bool selfGrows() { return false; }
	//virtual bool ignoresGrowthInhibition() { return false; }
	//virtual int getApopLiklihood() { return 0; }
	//virtual bool avoidsApoptosis() { return false; }
	//virtual bool ignoresTelomere() { return false; }
	//virtual bool deregulatesCellEnergetics() { return false; }
	//virtual int getImmuneKillLiklihood() { return 0; }
	//virtual bool avoidsImmunity() { return false; }
	//virtual void setAvoidApop(bool avoid);

	//Trying without inheritance...
	bool selfGrows() { return selfGrowth; }
	bool ignoresGrowthInhibition() { return ignoreGrowthInhibition; }
	int getApopLiklihood() { return apopLiklihood; }
	bool avoidsApoptosis() { return avoidApoptosis; }
	bool ignoresTelomere() { return ignoreTelomere; }
	bool angiogenesis() { return susAngio; }
	//bool deregulatesCellEnergetics() { return deregulateCellEnergetics; }
	int getImmuneDeathLiklihood() { return immuneDeathLiklihood; }
	bool avoidsImmunity() { return avoidImmunity; }
	bool genomeUnstable() { return genomicInstability; }

	//oxygen stuff
	double getRequiredOxygey() { return requiredOxy; }
	void setRequiredOxygen(double newOxy);
	int getNumMuts() { return numMuts; }
	void calcNumMuts();

	//Setters for hallmarks if cancerous
	void setSelfGrowth(bool self) { selfGrowth = self; }
	void setIgnoreGrowth(bool ignore) { ignoreGrowthInhibition = ignore; }
	void setApopLiklihood(int likely) { apopLiklihood = likely; }
	void setAvoidApop(bool avoid) { avoidApoptosis = avoid; }
	void setIgnoreTelo(bool ignoreT) { ignoreTelomere = ignoreT; }
	void setSustainedAngio(bool sus) { susAngio = sus; }
	//void setDeregulate(bool dereg) { deregulateCellEnergetics = dereg; }
	void setImmuneDeathLiklihood(int immuneDeathL) { immuneDeathLiklihood = immuneDeathL; }
	void setAvoidImmunity(bool avoidImmune) { avoidImmunity = avoidImmune; }
	void setGenomicInstability(bool genInstability) { genomicInstability = genInstability; }

private:
    //List of neighbours
    int telomere; //Telomore int.. initialized to 100
    int mutationRate; //chance of mutation (m in santos paper)
	int evadeApoptosis; //a cell with n hallmarks mutated has an extra n/e likelihood of dying each cell cycle (unless evades apoptosis is true)
	int geneticInstabilityFactor; //increase of base mutation rate by a fator of this for cells with geneticInstability
    int randomDeathRate; //chance of random death
	int competition; //chance of killing a neighbour

   //Neighbour array should be a list of pointers to cells..
	//maybe should make this stllist.. List[Cell *]
	Cell * neighbours[8];
	
	//Locations on the grid so we can point it out later
	int i;
	int j;
	int k; //For when we move to 3D

	int state;
    bool mutated; //true if cell becomes cancerous

	//Cancer cell properties... turned on if mutated becomes true
	bool selfGrowth; //Hallmark 1
	bool ignoreGrowthInhibition; //Hallmark 2
	int apopLiklihood; //cell with n hallmarks has an n/x liklihood of dying each cell cycle, unless evade apoptosis is on Hallmark 3
	bool avoidApoptosis;  //Hallmark 3
	bool ignoreTelomere; //Hallmark 4
	bool susAngio; //Hallmark 5
	bool avoidImmunity; //Emerging hallmark 2

	//Hallmark 6 is ignored because it is after avascular growth
	bool deregCellEnergetics; //Hallmark 7... not sure how to do this one
	int immuneDeathLiklihood; //cell with n hallmarks has n/x liklihook of being killed by immune system unless avoid immune is one //Hallmark 8
	bool genomicInstability; //Enabiling characteristic 2.
	int angiogenesisImmunity;
	int avoidImmuneParam;

	//Oxygen 
	bool consumed;
	double requiredOxy;
	int numMuts; //number of mutations a cell has 
	void calcRequiredOxy();
};
 
 
} /* namespace std */
#endif /* CELL_H_ */