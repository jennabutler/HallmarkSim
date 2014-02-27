/**
 *
 * \class BinaryFluid computes the flow of oxygen in blood
 * \description Standard Lattice Boltzmann Binary fluid a la Alex thesis
 * \note Compile line on Ultra sun: cc -fast -x05 binary.c -lm -o binary
 *
 * \note Compile line on DEC: cc -03 binary.c -lm -o binary
 */
/* Complie line on Ultra sun
cc -fast -xO5 binary.c -lm -o binary
 */
/* Compile line on DEC
cc -O3 binary.c -lm -o binary
 */
#include "stdafx.h"
#include "BinaryFluid.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

//#include <unistd.h>


//Commented these out
/*
void parametercalc(void);
void streamout(void);
void equilibriumdist(void);
void initialize(void);*/

using namespace std;


//Always correct this at the beginning
//#define L 784
#define L 900 //!< L is the size of the lattice grid
#define Nmax 10000

#define eachIteration 75 //!< eachIteration: is the number of times the fluid will update for each iteration of cancer growth

#define temperature 0.7
#define kappa 0.02
//#define diffusion_real 1.34e-4 //(mm^2/ms)
#define diffusion_real 2.0e-6 //(mm^2/ms) //Oct 16th 2012
#define lambda 1.1
#define viscosity 0.6e-9 //(kg/(ms*mm)

//NOTE:  The total density throughout the simulation domain is set to:
//       density = density_blood + 0.7*density_oxygen.
//       At the simulation boundary, the density difference is set equal to:
//       phi = density_blood - density_oxygen.
//       This gives a value for the density of oxygen at the boundary of:
//       0.85*density_oxygen.
#define density_oxygen_real 9.0e-9 //(mol/(mm^3))
#define density_blood_real 3.0e-7 //(mol/(mm^3))
#define tau2 4.0
#define dx 0.075
#define dt 3.0
#define dm 1.0e-10

double Gamma,diffusion;
double tau1;
double density_blood,density_oxygen;

#define Pi 3.141592653589793
#define TwoPi 6.283185307179586

double f1[L][L][9],g1[L][L][9],f2[L][L][9],g2[L][L][9],Ff[L][L][9];
double Fc[L][L][9],Gc[L][L][9];
double feq[L][L][9],geq[L][L][9];
double density[L][L],phi[L][L],u[L][L][2];
int e[9][2];
double (*f)[L][9],(*g)[L][9],(*fnew)[L][9],(*gnew)[L][9];

//My own varibales
int scaleFactor;
int oxygenGridSize;
int topBound;
int bottomBound;
int leftBound;
int rightBound;
int numberOfCells = 300; //2/3

/*
 * Default constructor
 */
BinaryFluid::BinaryFluid(){

	//Scale the input parameters to LB units.
	rescale();

	initialize();
	//printf("%f %f %i %i %i %i\n",1.0,1.0,L,L,1,2);
	//printf("vx\n");
	//printf("vy\n");
	//printf("phi\n");
}

/** \description Constructor for BinaryFluid
 * \param cancerGridSize: This size of one side of the cancer grid
 *
 * Constructor takes cancerGridSize to get the scale factor set up
 * Calls initialization of the fluid, and then runs one cycle of fluid update
 */
BinaryFluid::BinaryFluid(int cancerGridSize){
	scaleFactor = L/cancerGridSize;
	//scaleFactor=1; //2/3
	topBound = 0;
	bottomBound = L-1;
	leftBound = 0;
	rightBound = L-1;

	//Scale the input parameters to LB units
	rescale();

	//Initialize fluid
	initialize();
	//Update it to get Ff to 0
	updateBinaryFluid(); //2/10
	//Run it a few times to get oxygen into middle
	for (int i=0; i<10; i++){
		updateBinaryFluidSmall(); //2/10
	}
}

/**
 * \description Scales the input parameters to be in LB units
 */
void BinaryFluid::rescale(){
	//Scale the input parameters to LB units.
	density_oxygen = density_oxygen_real*dx*dx*dx/dm;
	density_blood = density_blood_real*dx*dx*dx/dm;

	diffusion = diffusion_real*dt/dx/dx;

	//I put in 1.05e-6 kg/mm^3 for the density of blood.
	tau1 = 3.0*viscosity*dt/1.05e-6/dx/dx + 0.5;

	Gamma = 2.0*density_blood*diffusion/(lambda*(tau2-0.5));
}


void BinaryFluid::updateBinaryFluidSmall(){
	int i,j,k,n;
	int iup,idwn,jup,jdwn;
	double (*tmp)[L][9];

	parametercalc();
	equilibriumdist();

	for (i=0; i<L; i++)
		for (j=0; j<L; j++)
			for (k=0; k<9; k++) {
				Fc[i][j][k]= (feq[i][j][k]-f[i][j][k])/tau1;
				Gc[i][j][k]= (geq[i][j][k]-g[i][j][k])/tau2;
			}

	for (i=0; i<L; i++) {
		if (i==L-1) iup=0; else iup=i+1;
		if (i==0) idwn=L-1; else idwn=i-1;
		for (j=0; j<L; j++) {
			if (j==L-1) jup=0; else jup=j+1;
			if (j==0) jdwn=L-1; else jdwn=j-1;

			fnew[i][j][2]=f[i][jdwn][2]+Fc[i][jdwn][2];
			fnew[i][j][0]=f[i][j][0]+Fc[i][j][0];
			fnew[i][j][4]=f[i][jup][4]+Fc[i][jup][4];
			fnew[i][j][5]=f[idwn][jdwn][5]+Fc[idwn][jdwn][5];
			fnew[i][j][1]=f[idwn][j][1]+Fc[idwn][j][1];
			fnew[i][j][8]=f[idwn][jup][8]+Fc[idwn][jup][8];
			fnew[i][j][6]=f[iup][jdwn][6]+Fc[iup][jdwn][6];
			fnew[i][j][3]=f[iup][j][3]+Fc[iup][j][3];
			fnew[i][j][7]=f[iup][jup][7]+Fc[iup][jup][7];

			gnew[i][j][2]=g[i][jdwn][2]+Gc[i][jdwn][2];
			gnew[i][j][0]=g[i][j][0]+Gc[i][j][0];
			gnew[i][j][4]=g[i][jup][4]+Gc[i][jup][4];
			gnew[i][j][5]=g[idwn][jdwn][5]+Gc[idwn][jdwn][5];
			gnew[i][j][1]=g[idwn][j][1]+Gc[idwn][j][1];
			gnew[i][j][8]=g[idwn][jup][8]+Gc[idwn][jup][8];
			gnew[i][j][6]=g[iup][jdwn][6]+Gc[iup][jdwn][6];
			gnew[i][j][3]=g[iup][j][3]+Gc[iup][j][3];
			gnew[i][j][7]=g[iup][jup][7]+Gc[iup][jup][7];
		}
	}
	tmp=f;
	f=fnew;
	fnew=tmp;
	tmp=g;
	g=gnew;
	gnew=tmp;
}


/**
 * \description This is the simulation step for the binary fluid; This is where it does the streaming and collision;
 * This is run "eachIteration" number of times for each iteration of cancer growth
 *
 */
void BinaryFluid::updateBinaryFluid(){

	int i,j,k,n;
	int iup,idwn,jup,jdwn;
	double (*tmp)[L][9];

	//initialize();
	//printf("%f %f %i %i %i %i\n",1.0,1.0,L,L,1,2);
	//printf("vx\n");
	//printf("vy\n");
	//printf("density\n");
	//printf("phi\n");


	for (n=0; n<eachIteration; n++) {

		for(i=0; i<L; i++)
			for(j=0; j<L; j++)
				for(k=0; k<9; k++)
					Ff[i][j][k] = 0.0;


		parametercalc();

		/*for (i=0; i<L; i++){
			for (j=0; j<L; j++){
				if (i>60 || i<100){
					if (j>60 || j<100){
						for(k=0; k<9; k++)
							Ff[i][j][k] = -0.0000005;
					}
				}
			}
		}*/



		equilibriumdist();

		for (i=0; i<L; i++)
			for (j=0; j<L; j++)
				for (k=0; k<9; k++) {
					Fc[i][j][k]= (feq[i][j][k]-f[i][j][k])/tau1;
					Gc[i][j][k]= (geq[i][j][k]-g[i][j][k])/tau2;
				}

		for (i=0; i<L; i++) {
			if (i==L-1) iup=0; else iup=i+1;
			if (i==0) idwn=L-1; else idwn=i-1;
			for (j=0; j<L; j++) {
				if (j==L-1) jup=0; else jup=j+1;
				if (j==0) jdwn=L-1; else jdwn=j-1;

				fnew[i][j][2]=f[i][jdwn][2]+Fc[i][jdwn][2];
				fnew[i][j][0]=f[i][j][0]+Fc[i][j][0];
				fnew[i][j][4]=f[i][jup][4]+Fc[i][jup][4];
				fnew[i][j][5]=f[idwn][jdwn][5]+Fc[idwn][jdwn][5];
				fnew[i][j][1]=f[idwn][j][1]+Fc[idwn][j][1];
				fnew[i][j][8]=f[idwn][jup][8]+Fc[idwn][jup][8];
				fnew[i][j][6]=f[iup][jdwn][6]+Fc[iup][jdwn][6];
				fnew[i][j][3]=f[iup][j][3]+Fc[iup][j][3];
				fnew[i][j][7]=f[iup][jup][7]+Fc[iup][jup][7];

				gnew[i][j][2]=g[i][jdwn][2]+Gc[i][jdwn][2];
				gnew[i][j][0]=g[i][j][0]+Gc[i][j][0];
				gnew[i][j][4]=g[i][jup][4]+Gc[i][jup][4];
				gnew[i][j][5]=g[idwn][jdwn][5]+Gc[idwn][jdwn][5];
				gnew[i][j][1]=g[idwn][j][1]+Gc[idwn][j][1];
				gnew[i][j][8]=g[idwn][jup][8]+Gc[idwn][jup][8];
				gnew[i][j][6]=g[iup][jdwn][6]+Gc[iup][jdwn][6];
				gnew[i][j][3]=g[iup][j][3]+Gc[iup][j][3];
				gnew[i][j][7]=g[iup][jup][7]+Gc[iup][jup][7];
			}
		}
		tmp=f;
		f=fnew;
		fnew=tmp;
		tmp=g;
		g=gnew;
		gnew=tmp;
	}

}

/**
 * \description Updates the equilibrium distribution for f and g
 */
void BinaryFluid::equilibriumdist(void)
{
	double A0,A1,A2,B1,B2,C0,C1,C2,D1,D2,G1xx,G2xx,G1xy,G2xy,G1yy,G2yy;
	double H0,H1,H2,K1,K2,J0,J1,J2,Q1,Q2;
	double dphidx,dphidy,rho,phiij,mu,usq,udote,laplacian;
	int i,j,k,iup,idwn,jup,jdwn;

	for (i=0; i<L; i++) {
		if (i==L-1) iup=0; else iup=i+1;
		if (i==0) idwn=L-1; else idwn=i-1;
		for (j=0; j<L; j++) {
			if (j==L-1) jup=0; else jup=j+1;
			if (j==0) jdwn=L-1; else jdwn=j-1;
			rho=density[i][j];
			phiij=phi[i][j];
			laplacian=kappa*((phi[iup][j]-2.0*phiij+phi[idwn][j])+
					(phi[i][jup]-2.0*phiij+phi[i][jdwn]));
			mu= -lambda/2.0*phiij/rho+
					temperature/2.0*log(((rho+phiij)/(rho-phiij)))-laplacian;
			A2= (rho*temperature-phiij*laplacian)/8.0;
			A1= 2.0*A2;
			A0= rho-12.0*A2;
			B2= rho/12.0;
			B1= 4.0*B2;
			C2= -rho/16.0;
			C1= 2.0*C2;
			C0= -3.0*rho/4.0;
			D2= rho/8.0;
			D1= 4.0*D2;
			dphidx=(phi[iup][j]-phi[idwn][j])/2.0;
			dphidy=(phi[i][jup]-phi[i][jdwn])/2.0;
			G2xx= kappa/16.0*(dphidx*dphidx-dphidy*dphidy);
			G2xy= kappa/8.0*dphidx*dphidy;
			G2yy= -G2xx;
			G1xx= 4.0*G2xx;
			G1xy= 4.0*G2xy;
			G1yy= 4.0*G2yy;
			H2= Gamma/8.0*mu;
			H1= 2.0*H2;
			H0= phiij-12.0*H2;
			K2= phiij/12.0;
			K1= 4.0*K2;
			J2= -phiij/16.0;
			J1= 2.0*J2;
			J0= -3.0*phiij/4.0;
			Q2= phiij/8.0;
			Q1= 4.0*Q2;

			usq=u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1];
			feq[i][j][0]=A0+C0*usq;
			geq[i][j][0]=H0+J0*usq+tau2*Ff[i][j][0];
			for (k=1; k<=4; k++) {
				udote=u[i][j][0]*e[k][0]+u[i][j][1]*e[k][1];
				feq[i][j][k]=A1+B1*udote+C1*usq+D1*udote*udote+G1xx*e[k][0]*e[k][0]+
						2.0*G1xy*e[k][0]*e[k][1]+G1yy*e[k][1]*e[k][1];
				geq[i][j][k]=H1+K1*udote+J1*usq+Q1*udote*udote+tau2*Ff[i][j][k];
			}
			for (k=5; k<=8; k++) {
				udote=u[i][j][0]*e[k][0]+u[i][j][1]*e[k][1];
				feq[i][j][k]=A2+B2*udote+C2*usq+D2*udote*udote+G2xx*e[k][0]*e[k][0]+
						2.0*G2xy*e[k][0]*e[k][1]+G2yy*e[k][1]*e[k][1];
				geq[i][j][k]=H2+K2*udote+J2*usq+Q2*udote*udote+tau2*Ff[i][j][k];
			}
		}
	}
}


void BinaryFluid::streamout(void)
{
	int i,j;

	/*   for(i=0; i<L; i++) { */
	/*     for (j=0; j<L; j++) { */
	/*       printf("%16.12f %16.12f %16.12f %16.12f\n", */
	/* 	     u[i][j][0],u[i][j][1],density[i][j],phi[i][j]); */
	/*     } */
	/*   } */

	cout << fixed;
	for(i=0; i<L; i++){
		for(j=0; j<L; j++){
			cout << setprecision(15) << phi[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

/**
 * \description calculates the parameters for the binary fluid; Sets the boundary conditions
 */
void BinaryFluid::parametercalc(void)
{
	int i,j,k;

	for (i=0; i<L; i++) {
		for (j=0; j<L; j++) {
			density[i][j]=0.0;
			phi[i][j]=0.0;
			u[i][j][0]=u[i][j][1]=0.0;
			for (k=0; k<9; k++) {
				density[i][j] += f[i][j][k];
				phi[i][j] += g[i][j][k];
				u[i][j][0] += f[i][j][k]*e[k][0];
				u[i][j][1] += f[i][j][k]*e[k][1];
			}
			//Have fixed boundary conditions
			if(i==0 || i == L-1 || j==0 || j==L-1)
				phi[i][j] = density_blood-density_oxygen;

			u[i][j][0]=u[i][j][0]/density[i][j];
			u[i][j][1]=u[i][j][1]/density[i][j];
		}
	}
}


void BinaryFluid::initialize(void)
{
	int i,j,k;

	e[0][0]= 0;
	e[0][1]= 0;
	e[1][0]= 1;
	e[1][1]= 0;
	e[2][0]= 0;
	e[2][1]= 1;
	e[3][0]= -1;
	e[3][1]= 0;
	e[4][0]= 0;
	e[4][1]= -1;
	e[5][0]= 1;
	e[5][1]= 1;
	e[6][0]= -1;
	e[6][1]= 1;
	e[7][0]= -1;
	e[7][1]= -1;
	e[8][0]= 1;
	e[8][1]= -1;

	f=f1;
	g=g1;
	fnew=f2;
	gnew=g2;
	for (i=0; i<L; i++) {
		for (j=0; j<L; j++) {
			//density[i][j]=density_blood+0.4*density_oxygen; oct 2 jcamer7
			density[i][j]=density_blood+0.7*density_oxygen;
			/*phi[i][j]= 0.0+0.04*(drand48()-0.5);*/
			u[i][j][0]=0.0;
			u[i][j][1]=0.0;
			/*    phi[i][j] = 1.12*(1.0+tanh((1.0*i-35.0)/4.0)-tanh((1.0*i-10.0)/4.0)); */
			/*       if (i== L/2 && j==L/2) phi[i][j]=0.1; */
			/*       if (i>L/2) phi[i][j]= 1.12 ; */
			/*if (i > 20 && i < 30 && j > 20 && j < 30){
				phi[i][j] = -0.05;
			}
			else{
				phi[i][j] = 0.05;
			}*/
			//Initialize boundary
			if(i==0 || i == L-1 || j==0 || j==L-1){
				phi[i][j] = density_blood-density_oxygen;
			}
			else {
				//phi[i][j] = density_blood-0.4*density_oxygen; oct 2nd jcamer7
				phi[i][j] = density_blood-0.7*density_oxygen;
			}


			/*       if ((i-L/2)*(i-L/2)+(j-L/2)*(j-L/2) < 64) */
			/* 	phi[i][j]=1.0; */
			/*       else */
			/* 	phi[i][j]= -1.0; */
			for (k=0; k<9; k++) {
				f[i][j][k]=density[i][j]/9.0;
				g[i][j][k]=phi[i][j]/9.0;
			}
		}
	}

	/*
  phi[L/2][L/2]=1.5;
  for (k=0; k<9; k++) {
    g[L/2][L/2][k]=phi[L/2][L/2]/9.0;
  }
	 */
}

/**
 * \description Method to get the oxygen values at a particular grid location
 * Must calculate the correct values based on a scale
 * Returns a sum of all the oxygen at a grid point on the oxygen grid for that point on the cancer grid
 */
double BinaryFluid::getOxygenValue(int i, int j){
	double oxygenValue = 0.0;

	//Loop through the oxygen grid points that make up a single cancer grid, and sum the oxygen at that point
	int oxyI = i*scaleFactor;
	int oxyJ = j*scaleFactor;
	int boundaryI = oxyI+scaleFactor;
	int boundaryJ = oxyJ + scaleFactor;
	for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
		for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
			oxygenValue = oxygenValue + 0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]); //updated Sept 26 2012
		}
	}
	return oxygenValue;
}


/**
 * \description Method to do a formatted print of the oxygen values
 */
string BinaryFluid::formatOxygenValue(){
	std::stringstream oxyValues;
	int boundaryI;
	int boundaryJ;
	int oxyI;
	int oxyJ;
	double oxyValue = 0;
	for (int i=0; i<numberOfCells; i++){
		for (int j=0; j<numberOfCells; j++){
			oxyI = i*scaleFactor;
			oxyJ = j*scaleFactor;
			boundaryI = oxyI+scaleFactor;
			boundaryJ = oxyJ + scaleFactor;
			oxyValue = 0;
			for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
				for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
					oxyValue = oxyValue + 0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]);
				}
			}
			oxyValues << oxyValue << " ";
		}
		oxyValues << endl;
	}
	return oxyValues.str();
}

/**
 * \description Method to take away oxygen at a cancer grid point
 * The i,j are the location of the CANCER cell so you need to really remove
 * the scaleFactor * i and j to get all the oxygen cells fixed up
 * \param amount this double must already be scaled per lattice boltzmann unit - not per cellular automaton unit. This occurs in Constants.h
 */
void BinaryFluid::consumeOxygen(int i, int j, double amount){
	//2/3
	//double currentOxygenAtLocation = 0.5*(density[i][j] - phi[i][j]);
	//bool enough = false;
	//if (currentOxygenAtLocation >= amount){
	//	enough = true;
	//}
	//if (enough){
	//	for (int k=0; k<9; k++){
	//		Ff[i][j][k] = amount/9.0; //nov 27
	//	}
	//}
	//else {
	//	//Don't remove the oxygen if not enough? Fix this //TO-DO
	//}
	//currentOxygenAtLocation = 0.5*(density[i][j] - phi[i][j]);

	
	int oxyI = i*scaleFactor;
	int oxyJ = j*scaleFactor;
	int boundaryI = oxyI + scaleFactor;
	int boundaryJ = oxyJ + scaleFactor;

	//oxygenAmount: array to hold the total amount of oxygen
	//              that needs to be removed from each LB grid site.
	double oxygenAmount[3][3];

	//hasOxy: set to: 0 if the site no longer has oxygen that can be removed.
	//              : 1 if oxygen is still available
	int hasOxy[3][3]; 

	double currentOxygenAtLocation;

	// First check to make sure that there is enough oxygen at each
	//  lattice-Boltzmann grid site.
	bool enough = true;
	for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
		for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
			currentOxygenAtLocation = 0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]);
			if (currentOxygenAtLocation < amount) {
				enough = false;
			}
		}
	}
	double oxygenTest1 =  0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]);
	
	// If there is enough oxygen at each site, evenly divide the total oxygen
	//  the cancer cell is consuming amoungst each LB site, and remove it:
	if(enough){
		for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
			for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
				for (int k=0; k<9; k++){
					double newAmount = amount/9.0;
					Ff[oxyI][oxyJ][k] = newAmount; //nov 27
				}
			}
		}
	}
	//else {
	//	cout << "not enough";
	//}
	//MARK not doing the weird tiny oxygen consumption thing for now

	// If there is not enough oxygen at each site, take away as much as possible
	//  from the sites that don't have enough, and then evenly distribute the 
	//  leftover amount amoungst the remaining sites (always checking to make sure a 
	//  site has enough oxygen available).
	else{

		//Set the total remaining oxygen that needs to be removed by the
		// cancer cell.
		double remainingOxy = amount*scaleFactor*scaleFactor;
		double remainingOxy_old;

		//Initialize the oxygenAmount array to 0.0, indicating that no oxygen has
		// been set aside to be removed yet.
		//Initialize the hasOxy array to 1, indicating that at each lattice-Boltzmann
		// grid point there is oxygen available to be removed.
		for(int m=0; m<scaleFactor; m++){
			for(int n=0; n<scaleFactor; n++){
				oxygenAmount[m][n] = 0.0;
				hasOxy[m][n] = 1;
			}
		}

		//Set the number of lattice Boltzmann sites that still have oxygen
		// that can be removed.
		int numOxy = scaleFactor*scaleFactor;
		int numOxy_old;

		while(remainingOxy > 1.0e-10){
			numOxy_old = numOxy;
			numOxy = 0;

			remainingOxy_old = remainingOxy;
			remainingOxy = 0.0;

			for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
				for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){

					//If the site still has available oxygen:
					if(hasOxy[oxyI-i*scaleFactor][oxyJ-j*scaleFactor] == 1){

						//Calculate the available oxygen (need to subtract off the
						// amount that already needs to be consumed).
						currentOxygenAtLocation = 0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]) -
								oxygenAmount[oxyI-i*scaleFactor][oxyJ-j*scaleFactor];

						//Want to remove the remaining amount of oxygen (remainingOxy_old),
						// evenly distributed amoungst the LB sites that still have oxygen available.

						//If there isn't enough oxygen at the site:
						if (currentOxygenAtLocation < remainingOxy_old/numOxy_old) {

							//Flag the site, remove all of the oxygen at the site, and calculate
							// the remaining oxygen that will now need to be removed in the next
							// iteration, by the remaining LB sites that still have oxygen.
							hasOxy[oxyI-i*scaleFactor][oxyJ-j*scaleFactor] = 0;
							oxygenAmount[oxyI-i*scaleFactor][oxyJ-j*scaleFactor] += currentOxygenAtLocation;
							remainingOxy += remainingOxy_old/numOxy_old - currentOxygenAtLocation;
						}

						//If there is enough oxygen at the site:
						else{

							//Remove the requested amount of oxygen, and add this site to the
							// count of LB sites that still have available oxygen.
							oxygenAmount[oxyI-i*scaleFactor][oxyJ-j*scaleFactor] += remainingOxy_old/numOxy_old;
							numOxy++;

						}
					}
				}
			}
		}

		// Remove the oxygen:
		for (oxyI=i*scaleFactor; oxyI < boundaryI; oxyI++){
			for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
				for(int k=0; k<9; k++){
					Ff[oxyI][oxyJ][k] = oxygenAmount[oxyI-i*scaleFactor][oxyJ-j*scaleFactor]/9.0;
				}
			}
		}

	}
}

/**
 * \description Method to output the oxygen values at each point of cancer cell grid
 * \note Here is where you specify where you want the files saved
 *
 */
std::string BinaryFluid::outputOxygenGrid(int iter){
	std::stringstream oxygenValues;
	double tempAmount = 0;
	for (int j=0; j<numberOfCells; j++){
		for (int i=0; i<numberOfCells; i++){
			int oxyI = i*scaleFactor;
			int oxyJ = j*scaleFactor;
			int boundaryI = oxyI + scaleFactor;
			int boundaryJ = oxyJ + scaleFactor;
			for (oxyI=i*scaleFactor; oxyI<boundaryI; oxyI++){
				for (oxyJ=j*scaleFactor; oxyJ<boundaryJ; oxyJ++){
					tempAmount = tempAmount + 0.5*(density[oxyI][oxyJ] - phi[oxyI][oxyJ]); //updated Sept 26 2012
				}
			}
			oxygenValues << ", " << setprecision(3) << tempAmount;
			tempAmount = 0;
		}
		oxygenValues << endl;
	}
//		for (int j=0; j<L; j++){
//			for (int i=0; i<L; i++){
//				tempAmount = 0.5*(density[i][j] - phi[i][j]);
//				oxygenValues << ", " << setprecision(3) << tempAmount;
//			}
//			oxygenValues << endl;
//		}

	std::string oxygenValuesString = oxygenValues.str();
	//char fileName[150];
//	sprintf(fileName, "/Users/jenna/Research/CancerModelAdditionalFiles/CAOutput/Dec2012/Dec17_2012Run1/OxygenFile%03d.txt", iter);
	//ofstream outfile;
	//outfile.open (fileName);
	//outfile << oxygenValuesString;
	//outfile.close();

	return oxygenValuesString;
}




/*
 * Print out oxygen
 */
void BinaryFluid::printOxygen(){
	for (int i=0; i<L; i++){
		for (int j=0; j<L; j++){
			cout << Ff[i][j][1] << " ";
		}
		cout << endl;
	}
}





