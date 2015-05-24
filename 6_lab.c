#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h>

float deltaX = 1;
float deltaY = 1;

float deltaTime = 0.2;
float time = 100;

size_t iterationAmount;

size_t sizeX = 128 * 4 - 1;
size_t sizeY = 128 * 4 - 1;

size_t nodesX;
size_t nodesY;

float defaultTemp;

double *mesh = NULL;

struct Edge {
    int type;
    float temp;
};

typedef struct Edge edge;

edge condition[4];

void initialize() {
    nodesX = sizeX / deltaX + 1;
    nodesY = sizeY / deltaY + 1;

    iterationAmount = time / deltaTime + 1;

    size_t i = 0;
    for (i = 0; i < 4; ++i) {
        condition[i].temp = 100;
        condition[i].type = 1;
    }
}

void createMesh(size_t myId) {
    if (!myId) {
        printf("I'm the main one\n");
        mesh = (double*) malloc(sizeof(double) * nodesX * nodesY);
    } else {
        mesh = NULL;
    }
}

void printSettings() {
    printf("sizeX = %d, sizeY = %d, deltaX = %f, deltaY = %f\n", (int)sizeX, (int)sizeY, deltaX, deltaY);
    printf("nodesX = %d, nodesY = %d\n", (int)nodesX, (int)nodesY);
    printf("Iteration amount = %d\n", (int)iterationAmount);
}

double cross(double top, double right, double bottom, double left, double central) {
    return central + deltaTime * (
            (left - 2 * central + right) / (deltaX * deltaX) +
            (top - 2 * central + bottom) / (deltaY * deltaY));
}
#define TESTPROC 0

void iterate(double *meshProc, size_t xx, double *left, double *right, double* result, int myId) {
    size_t i, j;
    for (j = 0; j < xx; ++j) {
        if (condition[2].type == 1) { // Первый род
            result[j * nodesY + 0] = condition[2].temp;
        } else { // Второй род
            result[j * nodesY + 0] = meshProc[j * nodesY + 1] - condition[2].temp * deltaY ;
        }
    }

    for (i = 1; i < nodesY - 1; ++i) {

        if (left == NULL) {
            if (condition[3].type == 1) { // Первый род
                result[0 * nodesY + i] = condition[3].temp;
            } else { // Второй род
                result[0 * nodesY + i] = meshProc[1 * nodesY + i] - condition[3].temp * deltaX ;
            }
        } else {
            j = 0;
            result[j * nodesY + i] = cross(meshProc[j * nodesY + i + 1], meshProc[(j + 1) * nodesY + i], meshProc[j * nodesY + i - 1], left[i], meshProc[j * nodesY + i]);
        }

        for (j = 1; j < xx - 1; ++j) {
            result[j * nodesY + i] = cross(meshProc[j * nodesY + i + 1], meshProc[(j + 1) * nodesY + i], meshProc[j * nodesY + i - 1], meshProc[(j - 1) * nodesY + i], meshProc[j* nodesY + i]);
        }

        j = xx - 1;
        if (right == NULL){
            if (condition[1].type == 1) { // Первый род
                result[(j) * nodesY + (i)] = condition[1].temp;
            } else { // Второй род
                result[(j) * nodesY + (i)] = meshProc[(j - 1) * nodesY + (i)] - condition[3].temp * deltaX ;
            }
        } else {
            result[(j) * nodesY + (i)] = cross(meshProc[(j) * nodesY + (i + 1)], right[i], meshProc[(j) * nodesY + (i - 1)], meshProc[(j - 1) * nodesY + (i)], meshProc[(j) * nodesY + (i)]);
        }


    }

    for (j = 0; j < xx; ++j) {
        if (condition[2].type == 1) { // Первый род
            result[(j) * nodesY + (nodesY - 1)] = condition[2].temp;
        } else { // Второй род
            result[(j) * nodesY + (nodesY - 1)] = meshProc[(j) * nodesY + (nodesY - 2)] - condition[2].temp * deltaY ;
        }
    }
}

void printTempField(double *toPrint) {
    size_t i, j;
    for ( i = 0; i < nodesY; ++i) {
        for (j = 0; j < nodesX; ++j) {
            printf("%5.1f ", toPrint[(j) * nodesY + (i)]);
        }
        putchar('\n');
    }
                    
}

void swap(double **a, double **b) {
    double *c = *a;
    *a = *b;
    *b = c;
}

int main(int argc, char* argv[]) {
    int myId;
    int procAmount;

    struct timeval tv1, tv2, dtv;
    struct timezone tz;

    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &procAmount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myId);

    initialize();
    size_t colPerProc = nodesX / procAmount;

    if (!myId) {
        printf("%d = totalProcAmount\n", procAmount);
        createMesh(myId);
        printSettings();
        printf("colPerProc = %d\n", (int)colPerProc);
    }


    double *pthreadCurrentData = (double*)malloc(sizeof(double) * colPerProc * nodesY);
    double *pthreadNextData = (double*)malloc(sizeof(double) * colPerProc * nodesY);
    double *left = NULL;
    if (myId != 0) {
        left = (double *) malloc(sizeof(double) * nodesY);
    }
    double *right = NULL;
    if (myId != procAmount - 1) {
        right = (double *) malloc(sizeof(double) * nodesY);
    }
    MPI_Scatter((void *)mesh, colPerProc * nodesY, MPI_DOUBLE, 
            (void *)pthreadCurrentData, colPerProc * nodesY, MPI_DOUBLE, 
            0 , MPI_COMM_WORLD);

    size_t iterationNum = 0;

    if (myId == 0) {
        gettimeofday(&tv1, &tz);
    }

    for (iterationNum = 0; iterationNum < iterationAmount; ++iterationNum ) {

        if ((myId % 2) == 1){ 
            if (myId > 0){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD);
            }
            if (myId + 1 < procAmount){
                MPI_Send((void *)&(pthreadCurrentData[(colPerProc - 1) * nodesY]), nodesY, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD);
            }   
            if (myId > 0){
                MPI_Recv((void *)left, nodesY, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD, &status);
            }
            if (myId + 1 < procAmount){
                MPI_Recv((void *)right, nodesY, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD, &status);
            }   
        }
        if (myId % 2 == 0){
            if (myId + 1 < procAmount){
                MPI_Recv((void *)right, nodesY, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD, &status);
            }

            if (myId > 0){
                MPI_Recv((void *)left, nodesY, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD, &status);
            }
            if (myId + 1 < procAmount){
                MPI_Send((void *)&(pthreadCurrentData[(colPerProc - 1) * nodesY]), nodesY, MPI_DOUBLE, myId + 1, 0, MPI_COMM_WORLD);
            }   
            if (myId > 0){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, 0, MPI_COMM_WORLD);
            }   
        }
        iterate(pthreadCurrentData, colPerProc, left, right, pthreadNextData, myId);
        if (iterationNum % 5 == 0) {
            //MPI_Gather((void*)pthreadNextData, colPerProc * nodesY, MPI_DOUBLE,
            //        (void*)mesh, colPerProc * nodesY, MPI_DOUBLE,
            //        0, MPI_COMM_WORLD);
            if (myId == 0 && 0) {
                printTempField(mesh);
                printf("IterationNum = %d\n", (int)iterationNum);
            }
        }
        swap(&pthreadCurrentData, &pthreadNextData);
    }
    if (myId == 0)
        gettimeofday(&tv2, &tz);
    free(pthreadNextData);
    free(pthreadCurrentData);
    if (!myId) {
        free(mesh);
    }
    if (myId == 0) {
        dtv.tv_sec = tv2.tv_sec - tv1.tv_sec;
        dtv.tv_usec = tv2.tv_usec - tv1.tv_usec;
        printf("Finish. Time to compute:\n");
        printf("%ld mc\n", dtv.tv_sec * 1000000 + dtv.tv_usec);
    }
    MPI_Finalize();
}
