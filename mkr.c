#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <sys/time.h>


int flag = 0;
float  a = 1;

float  deltaX = 1;
float  deltaY = 1;

double deltaTime = 0.0;
double  fullTime = 0.0;

int nodeCountX;
int nodeCountY;

int steps;
int perPthread;

double* top;
double* bottom;

double* nextData;
double* currentData;
double* pthreadNextData;
double* pthreadCurrentData;

struct Edge {
    int type;
    float temp;
};

int total;

void InitialFirsLayer(){
    int i, j, k;
    i = 0;
    for (j = 0; j < nodeCountY; ++j)
        for (k = 0; k < nodeCountX; ++k)
            currentData[j * nodeCountX + k] = ++i;
}


void Initial(int argc, char* argv[]){
    if (argc > 5) { 
        exit(0);
    }

    nodeCountY = atoi(argv[1]); // высота
    nodeCountX = atoi(argv[2]); // ширина
    fullTime = atof(argv[3]); // временной интервал
    deltaTime = atof(argv[4]);
    steps = (int)fullTime / deltaTime;

    currentData = (double*)malloc(nodeCountX * nodeCountY * sizeof(double));

    nextData = (double*)malloc(nodeCountX * nodeCountY * sizeof(double));

    pthreadCurrentData = (double*)malloc(nodeCountX * (perPthread + 2) * sizeof(double));

    pthreadNextData = (double*)malloc(nodeCountX * (perPthread + 2) * sizeof(double));

    bottom = (double*)malloc(nodeCountX * sizeof(double*));

    top = (double*)malloc(nodeCountX  * sizeof(double));
}



void GetEdge(struct Edge mas[]) {
    int i;
    FILE* fp = fopen("settings.txt", "r");
    float temp;
    for (i = 0; i < 4; ++i) {
        fscanf(fp, "%d %f\n", &(mas[i].type), &temp);
        mas[i].temp = (double) temp;
    }
    fclose(fp);   
}

//void SetEdge(double* pthreadCurrentData, double* pthreadNextData, int myrank){
void SetEdge(int myrank){
    int j, k;
    struct Edge edge[4];
    GetEdge(edge);
    //for (j = 0; j < 4; ++j)
    // printf("edge[%d] = %f\n", j, edge[j].temp);
    // printf("myrank = %d\n", myrank);
    if (edge[0].type == 1){ //left
        //printf("left\n");
        //printf("per = %d", perPthread);
        for (j = 0; j < perPthread; ++j){
            pthreadCurrentData[j * nodeCountX] = edge[0].temp;
            //printf ("%f = %f  ",  pthreadCurrentData[j * nodeCountX], edge[0].temp);
        }
    }
    else{
        for (j = 0; j < perPthread; ++j)
            pthreadCurrentData[j * nodeCountX] = pthreadCurrentData[j * nodeCountX + 1] + edge[0].temp * deltaX / a;
    }
    if (myrank == total - 1) {
        if (edge[1].type == 1){ // top
            // printf("top\n");
            for (k = 0; k < nodeCountX; ++k){
                pthreadCurrentData[(nodeCountY - 1) * nodeCountX + k] = edge[1].temp;
                //printf ("index = %d:  ", (nodeCountY - 1) * nodeCountX + k);
                //printf ("%f = %f\n",  pthreadCurrentData[(nodeCountY - 1) * nodeCountX + k], edge[1].temp);
            }
        }
        else{
            for (k = 0; k < nodeCountX; ++k)
                pthreadCurrentData[(nodeCountY - 1) * nodeCountX + k] = pthreadCurrentData[(nodeCountY - 2) * nodeCountX + k] + edge[1].temp * deltaY / a;
        }
    }

    if (edge[2].type == 1){ // right
        // printf("right\n");
        for (j = 0; j < perPthread; ++j)
            pthreadCurrentData[j * nodeCountX + nodeCountX - 1] = edge[2].temp;
    }
    else{
        for (j = 0; j < perPthread; ++j)
            pthreadCurrentData[j * nodeCountX + nodeCountX - 1] = pthreadCurrentData[j * nodeCountX + nodeCountX - 2]  + edge[2].temp * deltaX / a;
    }

    if (myrank == 0) {
        if (edge[3].type == 1){ //bottom
            // printf("bottom\n");
            for (k = 0; k < nodeCountX; ++k)
                pthreadCurrentData[k] = edge[3].temp;
        }
        else{
            for (k = 0; k < nodeCountX; ++k)
                pthreadCurrentData[k] = pthreadCurrentData[nodeCountX + k] + edge[3].temp * deltaY / a;
        }
    }
}


void Method(double* pthreadCurrentData, double* pthreadNextData, int rank){
    --rank;
    // Kostil
    
    int j, k;
    for(j = 1; j < nodeCountY - 1; ++j){
        for (k = 1; k < nodeCountX - 1; ++k){
            double tempX = a * (pthreadCurrentData[j * nodeCountX + k - 1] - 2 * pthreadCurrentData[j * nodeCountX + k] + pthreadCurrentData[j * nodeCountX + k + 1]) * deltaTime / (2 * deltaX * deltaX);
            double tempY = a * (pthreadCurrentData[(j - 1) * nodeCountX + k] - 2 * pthreadCurrentData[j * nodeCountX + k] + pthreadCurrentData[(j + 1) * nodeCountX + k]) * deltaTime / (2 * deltaY * deltaY);
            pthreadNextData[j * nodeCountX + k] = tempX + tempY + pthreadCurrentData[j * nodeCountX + k];  
        }
    }
    // SetEdge(pthreadCurrentData, pthreadNextData, rank);
}


void PrintMas(int st){
    int j, k;
    printf ("step = %d:\n", st);
    for (j = 0; j < nodeCountY; ++j){   
        for (k = 0; k < nodeCountX; ++k)
            printf("%f ", currentData[j * nodeCountX + k]); 
        printf ("\n");
    }
    printf ("\n\n");
}

void PrintToFile(int interationNumber) {
    char fileName[10];
    sprintf(fileName, "step%d", interationNumber);
    FILE* fp = fopen(fileName, "w");
    int i, j;
    for (i = 0; i < nodeCountY; ++i) {
        for (j = 0; j < nodeCountY; ++j) {
            fprintf(fp, "%d %d %f\n", j, i, nextData[i * nodeCountX + j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}


void PrintResult() {
    FILE* fp = fopen("result", "w");
    int i;
    fprintf(fp, "set pm3d\n");
    fprintf(fp, "set view map\n");
    for (i = 0; i < steps; ++i){
        fprintf(fp, "splot 'step%d' with pm3d\n", i);
        fprintf(fp, "pause %f\n", deltaTime);
    }
    fprintf(fp, "pause -1\n");
    fclose(fp);
}



int main(int argc, char* argv[]){

    if(argc < 5){
        printf("Usage: ./progName nodeCountX nodeCountY fullTime timeDelta\n");
        return -1;
    }

    int myrank;

    MPI_Status s1, s2;

    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &total);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

    Initial(argc, argv);
    int i, j, k;

    if (nodeCountY % total){
        printf("Mistake\n");
          MPI_Finalize();
        exit(-1);
    }

    perPthread = nodeCountY / total;
    //printf("per = %d myrank = %d\n", perPthread, myrank);
    if (!myrank){
        InitialFirsLayer();
    }

    for (i = 0; i < 1; ++i){
        //if (!myrank)
        //  printf("Lol step = %d:\n", i);
        MPI_Scatter((void *)currentData, perPthread * nodeCountX, MPI_DOUBLE, (void *)pthreadCurrentData, perPthread * nodeCountX,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);


        //if (!myrank)
        //  printf("Second step = %d:\n", i);
        if (myrank){ //так как нужны данные от соседей, то двумя след. функциями их получение и передача своих
            MPI_Sendrecv(&pthreadCurrentData[nodeCountX * (perPthread - 1)], nodeCountX, MPI_DOUBLE, myrank - 1, 0, bottom, nodeCountX, MPI_DOUBLE, myrank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &s1);
        }
        if(myrank != total - 1) {
            MPI_Sendrecv(pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank + 1, 1, top, nodeCountX, MPI_DOUBLE, myrank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &s2);
        }
        /*if (myrank % 2){ 
          if (myrank - 1 >= 0){
          MPI_Send((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank - 1, i, MPI_COMM_WORLD);
          }
          if (myrank + 1 < total){
          MPI_Send((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank + 1, i, MPI_COMM_WORLD);
          }   
          if (myrank - 1 >= 0){
          MPI_Recv((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank - 1, i, MPI_COMM_WORLD, &status);
          }
          if (myrank + 1 < total){
          MPI_Recv((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank + 1, i, MPI_COMM_WORLD, &status);
          }   
          }
          if (myrank % 2 == 0){
          if (myrank + 1 < total){
          MPI_Recv((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank + 1, i, MPI_COMM_WORLD, &status);
          }

          if (myrank - 1 >= 0){
          MPI_Recv((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank - 1, i, MPI_COMM_WORLD, &status);
          }

          if (myrank + 1 < total){
          MPI_Send((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank + 1, i, MPI_COMM_WORLD);
          }   
          if (myrank - 1 >= 0){
          MPI_Send((void *)pthreadCurrentData, nodeCountX, MPI_DOUBLE, myrank - 1, i, MPI_COMM_WORLD);
          }   
          }*/
        if (myrank == 1){
            printf("step = %d, myrank = %d:\n", i, myrank);
            //printf("X = %d, Y = %d\n", nodeCountX, nodeCountY);
            printf("top:\n");
            for (k = 0; k < nodeCountX; ++k)
                printf("%f ", top[k]);

            printf("\nbottom:\n");
            for (k = 0; k < nodeCountX; ++k)
                printf("%f ",bottom[k]);
            printf("\n\n");

            for (k = 0; k < perPthread; ++k){
                for (j = 0; j < nodeCountX; ++j){
                    //printf("d ");
                    printf("%f  ", pthreadCurrentData[k * (nodeCountX) + j]);
                }
                printf("\n");
            }
        }
        //Method(pthreadCurrentData, pthreadNextData, myrank);
        //PrintToFile(i);
        SetEdge(myrank);

        /*for(i=1; i<n; i++)//сбор всех резульатов в рутовый
          {
          MPI_Gather((void *)(res+i*m), m, MPI_DOUBLE, (void *)(Res+i*(X+1)), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
          }*/

        //printf("myrank = %d, i = %f\n", myrank, (pthreadCurrentData + nodeCountX));

        //MPI_Gather((void *)(pthreadCurrentData + nodeCountX), perPthread * nodeCountX, MPI_DOUBLE, (void *)currentData, perPthread * nodeCountX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather((void *)(pthreadCurrentData), perPthread * nodeCountX, MPI_DOUBLE, (void *)currentData, perPthread * nodeCountX, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //printf(" proc %d\n", myrank);
        if (!myrank){
            printf("step = %d:\n", i);
            //printf("X = %d, Y = %d\n", nodeCountX, nodeCountY);
            for (k = 0; k < nodeCountY; ++k){
                for (j = 0; j < nodeCountX; ++j){
                    //printf("d ");
                    printf("%f  ", currentData[k * (nodeCountX) + j]);
                }
                printf("\n");
            }
            //PrintToFile(i);

            // PrintMas(i);
            //currentData = nextData;
        }

    }
    if (!myrank) {
        free(currentData);
        free(nextData);
        free(pthreadCurrentData);
        free(pthreadNextData);
    }

    MPI_Finalize();
    //PrintResult();
    return 0;   
}
