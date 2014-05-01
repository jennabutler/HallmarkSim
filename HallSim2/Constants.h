#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const int RAND_SEED = 5678;
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

const int GROWTH_RADIUS = 40; //200
const int GROWTH_H = 150; //2525 
const int GROWTH_K = 150; //2525
const int GROWTH_RADIUS_2 = 45;
const int GROWTH_RADIUS_3 = 50;
const int GROWTH_RADIUS_4 = 55;
 
const int LEFT_GROWTH = 125; //2425
const int RIGHT_GROWTH = 200; //2600
const int TOP_GROWTH = 50; //1500
const int BOTTOM_GROWTH = GRID_SIDE;

//Oxygen consumption rates
const double OXYGEN_CON_HEALTHY = 0.05;
const double OXYGEN_CON_QUIS = OXYGEN_CON_HEALTHY/2;
const double OXYGEN_DEREG = 0;
const double OXYGEN_CON_CANCER_AGG1 = OXYGEN_CON_HEALTHY/1.1;
const double OXYGEN_CON_CANCER_AGG2 = OXYGEN_CON_HEALTHY/1.2;
const double OXYGEN_CON_CANCER_AGG3 = OXYGEN_CON_HEALTHY/1.3;

//April 28
//healhy = 0.05
//quis = 0.025
//dereg = 0.025
//agg1 = 0.045
//agg2 = 0.04
//agg3 = 0.035
const int MUTATION_RATE = 400;

const double END_PERCENT = 0.9999;
//0,0072



#endif

//On heap:
// allCells
// start
// firstEvent
// daughterCell