/* How to Run
Compile Using:
  gcc -Werror -Wall -O3 -lm -fopenmp n-body_std.c
Run Using: ./a.out [NumberOfInterations NumberOfBodies]
	./a.out 10000 200
	gnuplot plot3D.sh
For gprof:
	gcc -Werror -Wall -lm -pg n-body_std.c
	./a.out 10000 200
	gprof ./a.out > analysis.txt
	gprof ./a.out | ./gprof2dot.py | dot -Tpng -o gprof_output.png
For perf:
	 perf record -g -- ./a.out
	 perf script | c++filt | ./gprof2dot.py -f perf | dot -Tpng -o perf_output.png

Code Ref:https://rosettacode.org/wiki/N-body_problem#C
*/
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "support.h"
// Chnage this value to reflect your ID number
#define ID 1011182
typedef struct
{
	double x, y, z;
} vector;

int bodies, timeSteps;
int SimulationTime = 0;
double *masses, GravConstant;
vector *positions, *velocities, *accelerations;

double getRND()
{
	return (100.0 - rand() % 200) / 50.0;
}
void initiateSystemRND(int bodies)
{
	int i;
	srand(ID);
	masses = (double *)malloc(bodies * sizeof(double));
	positions = (vector *)malloc(bodies * sizeof(vector));
	velocities = (vector *)malloc(bodies * sizeof(vector));
	accelerations = (vector *)malloc(bodies * sizeof(vector));
	GravConstant = 0.01;
	for (i = 0; i < bodies; i++)
	{
		masses[i] = 0.4; //(rand()%100)/100.0;//0.4 Good Example
		positions[i].x = getRND();
		positions[i].y = getRND();
		positions[i].z = getRND();
		velocities[i].x = getRND() / 5.0; // 0.0;
		velocities[i].y = getRND() / 5.0; // 0.0;
		velocities[i].z = getRND() / 5.0; // 0.0;
	}
}

void bubblesort(int t1[], int t2[], int n)
{
	int i, j, temp;
	for (i = 0; i < n - 1; i++)
	{
		for (j = 0; j < n - i - 1; j++)
		{
			if (t1[j] > t1[j + 1])
			{
				// swap elements in t1
				temp = t1[j];
				t1[j] = t1[j + 1];
				t1[j + 1] = temp;

				// swap corresponding elements in t2
				temp = t2[j];
				t2[j] = t2[j + 1];
				t2[j + 1] = temp;
			}
			else if (t1[j] == t1[j + 1])
			{
				if (t2[j] > t2[j + 1])
				{
					// swap elements in t1
					temp = t1[j];
					t1[j] = t1[j + 1];
					t1[j + 1] = temp;

					// swap corresponding elements in t2
					temp = t2[j];
					t2[j] = t2[j + 1];
					t2[j + 1] = temp;
				}
			}
		}
	}
}

void resolveCollisions()
{
	int i, j;
	double dx, dy, dz, md;

	for (i = 0; i < bodies - 1; i++)
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions[i].x - positions[j].x);
			dy = fabs(positions[i].y - positions[j].y);
			dz = fabs(positions[i].z - positions[j].z);
			// if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){
			if (dx < md && dy < md && dz < md)
			{
				vector temp = velocities[i];
				velocities[i] = velocities[j];
				velocities[j] = temp;
			}
		}
}


void resolveCollisions_static(int number_of_threads)
{
	int *t1 = (int *)malloc(50000 * sizeof(int));
	int *t2 = (int *)malloc(50000 * sizeof(int));

	int i, j, cnt = 0;
	double dx, dy, dz, md;
#pragma omp parallel for private(j, dx, dy, dz, md) schedule(static, bodies / number_of_threads) num_threads(number_of_threads)
	for (i = 0; i < bodies - 1; i++)
	{
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions[i].x - positions[j].x);
			dy = fabs(positions[i].y - positions[j].y);
			dz = fabs(positions[i].z - positions[j].z);
			// if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){

			if (dx < md && dy < md && dz < md)
			{
#pragma omp critical
				{

					*(t1 + cnt) = i;
					*(t2 + cnt) = j;
					cnt++;
				}
			}
		}
	}
	// printf(" %d",cnt);
	bubblesort(t1, t2, cnt);
	for (i = 0; i < cnt; i++)
	{
		// Swap Velocities
		vector temp = velocities[t1[i]];
		velocities[t1[i]] = velocities[t2[i]];
		velocities[t2[i]] = temp;
	}
}

void resolveCollisions_dynamic(int number_of_threads)
{
	int *t1 = (int *)malloc(50000 * sizeof(int));
	int *t2 = (int *)malloc(50000 * sizeof(int));

	int i, j, cnt = 0;
	double dx, dy, dz, md;
#pragma omp parallel for private(j, dx, dy, dz, md) schedule(dynamic, bodies / number_of_threads) num_threads(number_of_threads)
	for (i = 0; i < bodies - 1; i++)
	{
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions[i].x - positions[j].x);
			dy = fabs(positions[i].y - positions[j].y);
			dz = fabs(positions[i].z - positions[j].z);
			// if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){

			if (dx < md && dy < md && dz < md)
			{
#pragma omp critical
				{

					*(t1 + cnt) = i;
					*(t2 + cnt) = j;
					cnt++;
				}
			}
		}
	}
	// printf(" %d",cnt);
	bubblesort(t1, t2, cnt);
	for (i = 0; i < cnt; i++)
	{
		// Swap Velocities
		vector temp = velocities[t1[i]];
		velocities[t1[i]] = velocities[t2[i]];
		velocities[t2[i]] = temp;
	}
}

void resolveCollisions_guided(int number_of_threads)
{
	int *t1 = (int *)malloc(50000 * sizeof(int));
	int *t2 = (int *)malloc(50000 * sizeof(int));

	int i, j, cnt = 0;
	double dx, dy, dz, md;
#pragma omp parallel for private(j, dx, dy, dz, md) schedule(guided, bodies / number_of_threads) num_threads(number_of_threads)
	for (i = 0; i < bodies - 1; i++)
	{
		for (j = i + 1; j < bodies; j++)
		{
			md = masses[i] + masses[j];
			dx = fabs(positions[i].x - positions[j].x);
			dy = fabs(positions[i].y - positions[j].y);
			dz = fabs(positions[i].z - positions[j].z);
			// if(positions[i].x==positions[j].x && positions[i].y==positions[j].y && positions[i].z==positions[j].z){

			if (dx < md && dy < md && dz < md)
			{
#pragma omp critical
				{

					*(t1 + cnt) = i;
					*(t2 + cnt) = j;
					cnt++;
				}
			}
		}
	}
	// printf(" %d",cnt);
	bubblesort(t1, t2, cnt);
	for (i = 0; i < cnt; i++)
	{
		// Swap Velocities
		vector temp = velocities[t1[i]];
		velocities[t1[i]] = velocities[t2[i]];
		velocities[t2[i]] = temp;
	}
}

void computeAccelerations()
{
	int i, j;
	for (i = 0; i < bodies; i++)
	{
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (j = 0; j < bodies; j++)
		{
			if (i != j)
			{
				// accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				vector sji = {positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sji.x, s * sji.y, s * sji.z};
				accelerations[i].x += S.x;
				accelerations[i].y += S.y;
				accelerations[i].z += S.z;
			}
		}
	}
}

void computeAccelerations_static(int thread_num)
{
	int i, j;
#pragma omp parallel for private(i, j) shared(accelerations) schedule(static, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (j = 0; j < bodies; j++)
		{
			if (i != j)
			{
				// accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				vector sji = {positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sji.x, s * sji.y, s * sji.z};
				accelerations[i].x += S.x;
				accelerations[i].y += S.y;
				accelerations[i].z += S.z;
			}
		}
	}
}

void computeAccelerations_dynamic(int thread_num)
{
	int i, j;
#pragma omp parallel for private(i, j) shared(accelerations) schedule(dynamic, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (j = 0; j < bodies; j++)
		{
			if (i != j)
			{
				// accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				vector sji = {positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sji.x, s * sji.y, s * sji.z};
				accelerations[i].x += S.x;
				accelerations[i].y += S.y;
				accelerations[i].z += S.z;
			}
		}
	}
}

void computeAccelerations_guided(int thread_num)
{
	int i, j;
#pragma omp parallel for private(i, j) shared(accelerations) schedule(guided, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (j = 0; j < bodies; j++)
		{
			if (i != j)
			{
				// accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(mod(subtractVectors(positions[i],positions[j])),3),subtractVectors(positions[j],positions[i])));
				vector sij = {positions[i].x - positions[j].x, positions[i].y - positions[j].y, positions[i].z - positions[j].z};
				vector sji = {positions[j].x - positions[i].x, positions[j].y - positions[i].y, positions[j].z - positions[i].z};
				double mod = sqrt(sij.x * sij.x + sij.y * sij.y + sij.z * sij.z);
				double mod3 = mod * mod * mod;
				double s = GravConstant * masses[j] / mod3;
				vector S = {s * sji.x, s * sji.y, s * sji.z};
				accelerations[i].x += S.x;
				accelerations[i].y += S.y;
				accelerations[i].z += S.z;
			}
		}
	}
}

void computeVelocities()
{
	int i;
	for (i = 0; i < bodies; i++)
	{
		// velocities[i] = addVectors(velocities[i],accelerations[i]);
		vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
		velocities[i] = ac;
	}
}

void computeVelocities_static(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(velocities) schedule(static, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// velocities[i] = addVectors(velocities[i],accelerations[i]);
		vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
		velocities[i] = ac;
	}
}

void computeVelocities_dynamic(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(velocities) schedule(dynamic, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// velocities[i] = addVectors(velocities[i],accelerations[i]);
		vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
		velocities[i] = ac;
	}
}

void computeVelocities_guided(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(velocities) schedule(guided, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// velocities[i] = addVectors(velocities[i],accelerations[i]);
		vector ac = {velocities[i].x + accelerations[i].x, velocities[i].y + accelerations[i].y, velocities[i].z + accelerations[i].z};
		velocities[i] = ac;
	}
}

void computePositions()
{
	int i;
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}

void computePositions_static(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(positions) schedule(static, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}

void computePositions_dynamic(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(positions) schedule(dynamic, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}

void computePositions_guided(int thread_num)
{
	int i;
#pragma omp parallel for private(i) shared(positions) schedule(guided, bodies / thread_num) num_threads(thread_num)
	for (i = 0; i < bodies; i++)
	{
		// positions[i] = addVectors(positions[i],addVectors(velocities[i],scaleVector(0.5,accelerations[i])));
		vector sc = {0.5 * accelerations[i].x, 0.5 * accelerations[i].y, 0.5 * accelerations[i].z};
		vector ac = {velocities[i].x + sc.x, velocities[i].y + sc.y, velocities[i].z + sc.z};
		vector bc = {positions[i].x + ac.x, positions[i].y + ac.y, positions[i].z + ac.z};
		positions[i] = bc;
	}
}

void simulate()
{
	SimulationTime++;
	computeAccelerations();
	computePositions();
	computeVelocities();
	resolveCollisions();
}

void simulate_static(int threads)
{
	// I dont run with the implemented computePositions and computeVelocities because the overheads outweight the parallel benefits
	SimulationTime++;
	computeAccelerations_static(threads);
	computePositions();	 // computePositions_static(int threads)
	computeVelocities(); // computeVelocities_static(int threads)
	resolveCollisions_static(threads);
}

void simulate_dynamic(int threads)
{
	// I dont run with the implemented computePositions and computeVelocities because the overheads outweight the parallel benefits
	SimulationTime++;
	computeAccelerations_dynamic(threads);
	computePositions();	 // computePositions_dynamic(int threads)
	computeVelocities(); // computeVelocities_dynamic(int threads)
	resolveCollisions_dynamic(threads);
}

void simulate_guided(int threads)
{
	// I dont run with the implemented computePositions and computeVelocities because the overheads outweight the parallel benefits
	SimulationTime++;
	computeAccelerations_guided(threads);
	computePositions();	 // computePositions_guided(int threads)
	computeVelocities(); // computeVelocities_guided(int threads)
	resolveCollisions_guided(threads);
}

void printBodiesInfo(FILE *lfp, FILE *dfp)
{
	int j;
	for (j = bodies - 10; j < bodies; j++)
		fprintf(lfp, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions[j].x, positions[j].y, positions[j].z, velocities[j].x, velocities[j].y, velocities[j].z);
	fprintf(lfp, "-------------------------------------------------------------------------------------------\n");
	for (j = bodies - 10; j < bodies; j++)
		fprintf(stdout, "Body%d %f\t: %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j + 1, masses[j], positions[j].x, positions[j].y, positions[j].z, velocities[j].x, velocities[j].y, velocities[j].z);
	fprintf(stdout, "-------------------------------------------------------------------------------------------\n");
}

int main(int argc, char *argv[])
{
	int i;
	int threads = 1;
	FILE *lfp = fopen("./outputs/logfile.txt", "w");
	FILE *dfp = fopen("./outputs/data.dat", "w");
	if (lfp == NULL || dfp == NULL)
	{
		printf("Please create the ./outputs directory\n");
		return -1;
	}
	if (argc == 4)
	{
		timeSteps = atoi(argv[1]);
		bodies = atoi(argv[2]);
		threads = atoi(argv[3]);
	}
	else
	{
		printf("%%*** RUNNING WITH DEFAULT VALUES ***\n");
		timeSteps = 10000;
		bodies = 200;
	}
	initiateSystemRND(bodies);
	// fprintf(stdout,"Running With %d Bodies for %d timeSteps. Initial state:\n",bodies,timeSteps);
	// fprintf(stderr,"Running With %d Bodies for %d timeSteps. Initial state:\n",bodies,timeSteps);
	// fprintf(lfp,"Running With %d Bodies for %d timeSteps. Initial state:\n",bodies,timeSteps);
	// fprintf(lfp,"Body   \t\t\t:\t\tx\t\ty\t\t\tz\t\t|\t\tvx\t\t\tvy\t\t\tvz\t\t\n");
	// initiateSystem("input.txt");

	////////////////////////////////////////////////////////////////////Initial
	 printBodiesInfo(lfp, dfp);
	printf("Threads: %d\n", threads);
	////////////////////////////////////////////////////////////////////Serial
	// initiateSystemRND(bodies);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		simulate();
	}
	stopTime(0);

	printf("\nSimulation Time:");
	elapsedTime(0);
	 printBodiesInfo(lfp, dfp);
	////////////////////////////////////////////////////////////////////Static
	initiateSystemRND(bodies);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		simulate_static(threads);
	}
	stopTime(0);
	printf("\nSimulation Time:");
	elapsedTime(0);
	 printBodiesInfo(lfp, dfp);
	////////////////////////////////////////////////////////////////////Dynamic
	initiateSystemRND(bodies);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		simulate_dynamic(threads);
	}
	stopTime(0);
	printf("\nSimulation Time:");
	elapsedTime(0);
	 printBodiesInfo(lfp, dfp);
	////////////////////////////////////////////////////////////////////Guided
	initiateSystemRND(bodies);
	startTime(0);
	for (i = 0; i < timeSteps; i++)
	{
		simulate_guided(threads);
	}
	stopTime(0);
	printf("\nSimulation Time:");
	elapsedTime(0);
	 printBodiesInfo(lfp, dfp);
	////////////////////////////////////////////////////////////////////
	fclose(lfp);
	fclose(dfp);
	return 0;
}
