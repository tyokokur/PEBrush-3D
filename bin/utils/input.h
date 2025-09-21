/* input.h */
// Parameters that change from code to code //

// Change these
double wand   = 1e-01;
double wcmp   = 1e-01;
int and_NrMax =     3;
int maxIter    =  500;
double wopt = 1e-02;

// Box size of the seed
int win_nx = 25;
int win_ny = 25;
int win_nz = 250;

// Seed strategy
int init_opt = 1;

// System
double lx = 5.00;
double ly = 5.00;
double lz = 50.0;
double DGC_lzmax = 45.0;

// Files
char seq_file[30] = "../ac.txt";
char Win[30] = "../W_in.dat";


/*****************************
    IGNORE EVERYTHING BELOW 
*****************************/

// Solution
double c0_plus  = 10e-03;
double c0_minus = 10e-03;

// Block Polymer 
double b0 = 1.00;
double v0 = 4.19; // Add decimal points

// Initialization
double xCmax   = 0.60;
double xCneck  = 10;
double xCltot  = 0; 

// Convergence parameters
int and_it    =   1;
int and_rep_N =  -1;   
double wopt_PB = 1e-02;
double Sm_PB   = 1e-10;
int MaxIT_PB   = 800;

double Sm1 = 1e-09;
double Sm2 = 1e-05;
double Sm3 = 1e-05;

int HALF_BOOL  = 1;
int win_box    = 1;
