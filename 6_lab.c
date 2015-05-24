#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h>

float deltaX = 1;
float deltaY = 1;

float deltaTime = 0.2;
float time = 10;

size_t iterationAmount;

size_t sizeX = 7;
size_t sizeY = 7;

size_t nodesX;
size_t nodesY;

float defaultTemp;

double **mesh = NULL;
double **result = NULL;

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
        mesh = (double**) malloc(sizeof(double*) * nodesX);
        result = (double**) malloc(sizeof(double*) * nodesX);
        size_t i = 0;
        for (; i < nodesX; ++i) {
            mesh[i] = (double*) malloc (sizeof(double) * nodesY);
            result[i] = (double*) malloc (sizeof(double) * nodesY);
            size_t j;
            for (j = 0; j < nodesY; ++j) {
                mesh[i][j] = defaultTemp;
                result[i][j] = defaultTemp;
            }
        }
    } else {
        mesh = NULL;
        result = NULL;
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

void iterate(double *meshProc, size_t xx, double *left, double *right, double* result) {
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
                result[j][i] = condition[1].temp;
            } else { // Второй род
                result[j][i] = meshProc[j - 1][i] - condition[3].temp * deltaX ;
            }
        } else {
            result[j][i] = cross(meshProc[j][i + 1], right[i], meshProc[j][i - 1], meshProc[j - 1][i], meshProc[j][i]);
        }
    }

    for (j = 0; j < xx; ++j) {
        if (condition[2].type == 1) { // Первый род
            result[j][nodesY - 1] = condition[2].temp;
        } else { // Второй род
            result[j][nodesY - 1] = meshProc[j][nodesY - 2] - condition[2].temp * deltaY ;
        }
    }
}

void printTempField(double **toPrint) {
    size_t i, j;
    for ( i = 0; i < nodesY; ++i) {
        for (j = 0; j < nodesX; ++j) {
            printf("%5.1f ", toPrint[j][i]);
        }
        putchar('\n');
    }
                    
}

void swap(double ***a, double ***b) {
    double **c = *a;
    *a = *b;
    *b = c;
}

void MPI_Init(int *, char***);
int MPI_COMM_WORLD;
int MPI_DOUBLE;
typedef int MPI_Status;
void MPI_Comm_size(int, size_t*);
void MPI_Comm_rank(int, size_t*);

void MPI_Send(void *, size_t, int, size_t, size_t, int);
void MPI_Recv(void *, size_t, int, size_t, size_t, int, MPI_Status);

int main(int argc, char* argv[]) {
    size_t myId;
    size_t procAmount;
    MPI_Status s1, s2, status;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &procAmount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myId);

    if (!myId) {
        initialize();
        createMesh(myId);
        printTempField(mesh);
    }

    if (!myId) {
        printSettings();
    }
    size_t iterationNum = 0;

    size_t colPerProc = nodesX / procAmount;
    printf("colPerProc = %d\n", (int)colPerProc);
    double *pthreadCurrentData = (double*)malloc(sizeof(double) * nodesX * nodesY);

    for (iterationNum = 0; iterationNum < iterationAmount; ++iterationNum ) {

        if (myId % 2){ 
            if (myId > 0){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, i, MPI_COMM_WORLD);
            }
            if (myId + 1 < procAmount){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId + 1, i, MPI_COMM_WORLD);
            }   
            if (myId > 0){
                MPI_Recv((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, i, MPI_COMM_WORLD, &status);
            }
            if (myId + 1 < procAmount){
                MPI_Recv((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId + 1, i, MPI_COMM_WORLD, &status);
            }   
        }
        if (myId % 2 == 0){
            if (myId + 1 < procAmount){
                MPI_Recv((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId + 1, i, MPI_COMM_WORLD, &status);
            }

            if (myId > 0){
                MPI_Recv((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, i, MPI_COMM_WORLD, &status);
            }

            if (myId + 1 < procAmount){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId + 1, i, MPI_COMM_WORLD);
            }   
            if (myId > 1){
                MPI_Send((void *)pthreadCurrentData, nodesY, MPI_DOUBLE, myId - 1, i, MPI_COMM_WORLD);
            }   
        }

/*
        double **secondPart = &(mesh[1 * colPerProc]);
        double **secondPartResult = &(result[1 * colPerProc]);
        double **thirdPart= &(mesh[2 * colPerProc]);
        double **thirdPartResult = &(result[2 * colPerProc]);
        double **fourthPart = &(mesh[3 * colPerProc]);
        double **fourthPartResult = &(result[3 * colPerProc]);

        iterate(mesh, colPerProc, NULL, secondPart[0], result);
        iterate(secondPart, colPerProc, mesh[colPerProc - 1], thirdPart[0], secondPartResult);
        iterate(thirdPart, colPerProc, secondPart[colPerProc - 1], fourthPart[0], thirdPartResult);
        iterate(fourthPart, colPerProc, thirdPart[colPerProc - 1], NULL, fourthPartResult);
        putchar('\n');
        printTempField(result);
        swap(&mesh, &result);*/
    }

}
