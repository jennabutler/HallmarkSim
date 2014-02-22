#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const int RAND_SEED = 1234;
const int NUM_NEIGHBOURS = 8;
const int GRID_SIDE = 300; //5001
const int DIM = 2;
const int RUN = 1;

const int LEFT_BLOOD = 100; //2450 
const int RIGHT_BLOOD = 200; //2550
const int TOP_BLOOD = 0;
const int BOTTOM_BLOOD = 200; //2700

const int LEFT_BLOOD_2 = 75;
const int RIGHT_BLOOD_2 = 225;
const int TOP_BLOOD_2 = 0;
const int BOTTOM_BLOOD_2 = 225;

const int LEFT_BLOOD_3 = 50;
const int RIGHT_BLOOD_3 = 250;
const int TOP_BLOOD_3 = 0;
const int BOTTOM_BLOOD_3 = 250;

const int BLOOD_RADIUS = 40; //300
const int BLOOD_RADIUS_2 = 45;
const int BLOOD_RADIUS_3 = 50;
const int BLOOD_RADIUS_4 = 55;
const int BLOOD_H = 150; //2480
const int BLOOD_K = 150; //2480

const int GROWTH_RADIUS = 75; //200
const int GROWTH_H = 150; //2525 
const int GROWTH_K = 150; //2525
const int GROWTH_RAIDUS_2 = 80;
const int GROWTH_RADIUS_3 = 85;
const int GROWTH_RADIUS_4 = 90;
 
const int LEFT_GROWTH = 125; //2425
const int RIGHT_GROWTH = 200; //2600
const int TOP_GROWTH = 50; //1500
const int BOTTOM_GROWTH = GRID_SIDE;

const double OXYGEN_CONSUMPTION = 0.011;

const int MUTATION_RATE = 500;

const double END_PERCENT = 0.99;
//0,0072



#endif

//On heap:
// allCells
// start
// firstEvent
// daughterCell