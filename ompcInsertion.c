#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <omp.h>
#include <float.h>
#include <time.h>

#ifndef coordReader_H_
#define coordReader_H_

int readNumOfCoords();
double **readCoords();
void *writeTourToFile();
#endif
// Provided functions
int readNumOfCoords(char *filename);
double **readCoords(char *filename, int numOfCoords);
void *writeTourToFile(int *tour, int tourLength, char *filename);




// Function to calculate the Euclidean distance between two points
double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

// Function to generate a distance matrix from coordinates



double **generateDistanceMatrix(double **coords, int numOfCoords) {
    double **matrix = (double **)malloc(numOfCoords * sizeof(double *));
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i] = (double *)malloc(numOfCoords * sizeof(double));
    }
    #pragma omp parallel for
    for (int i = 0; i < numOfCoords; i++) {
        matrix[i][i] = 0; // Distance from a point to itself is 0
        for (int j = i + 1; j < numOfCoords; j++) {
            matrix[i][j] = distance(coords[i][0], coords[i][1], coords[j][0], coords[j][1]);
            matrix[j][i] = matrix[i][j]; // Use symmetry, avoid redundant calculation
        }
    }
    return matrix;
}

    void findNearestNeighbor(double **distanceMatrix, int numOfCoords, int currentNode, int *nearestIndex, double *nearestDistance) {
    double globalNearestDistance = DBL_MAX;
    int globalNearestIndex = -1;

    #pragma omp parallel
    {
        double localNearestDistance = DBL_MAX;
        int localNearestIndex = -1;

        #pragma omp for
        for (int i = 0; i < numOfCoords; i++) {
            if (i != currentNode) {
                double dist = distanceMatrix[currentNode][i];
                if (dist < localNearestDistance) {
                    localNearestDistance = dist;
                    localNearestIndex = i;
                }
            }
        }

        #pragma omp critical
        {
            if (localNearestDistance < globalNearestDistance) {
                globalNearestDistance = localNearestDistance;
                globalNearestIndex = localNearestIndex;
            }
        }
    }

    *nearestDistance = globalNearestDistance;
    *nearestIndex = globalNearestIndex;
}
 
int *cheapestInsertion(double **distanceMatrix, int numOfCoords) {
    int *tour = malloc((numOfCoords + 1) * sizeof(int));
    int *unvisited = malloc(numOfCoords * sizeof(int));
    // Initialize unvisited nodes
    for (int i = 1; i < numOfCoords; i++) {
        unvisited[i - 1] = i;
    }
    int unvisitedCount = numOfCoords - 1;

    // Start with vertex 0
    tour[0] = 0;
    int tourSize = 1; 

    // Find the nearest neighbor to 0 and add to the tour
    double nearestDistance;
    int nearestIndex;
    findNearestNeighbor(distanceMatrix, numOfCoords, 0, &nearestIndex, &nearestDistance); 
   
    tour[1] = nearestIndex; 
    tour[2] = 0; 
    tourSize = 3; 

    // Remove the nearest neighbor from the unvisited list
    for (int i = 0; i < unvisitedCount; i++) {
        if (unvisited[i] == nearestIndex) {
            unvisited[i] = unvisited[unvisitedCount - 1];
            unvisitedCount--;
            break;
        }
    }
    double globalminCost;
    int globalminCostIndex, globalinsertPosition;

    int numThreads = omp_get_max_threads();
    double *localMinCosts = malloc(numThreads * sizeof(double));
    int *localMinCostIndex = malloc(numThreads * sizeof(int));
    int *localInsertPositions = malloc(numThreads * sizeof(int));

    while (tourSize < numOfCoords + 1) {
        globalminCost = DBL_MAX;
        globalminCostIndex = -1;
        globalinsertPosition = -1;

        #pragma omp parallel
        {
            /* double localminCost = DBL_MAX;
            int localminCostIndex = -1, localinsertPosition = -1; */

            int threadID = omp_get_thread_num();
            localMinCosts[threadID] = DBL_MAX;
            localMinCostIndex[threadID] = -1;
            localInsertPositions[threadID] = -1;

        #pragma omp for
        for (int idx = 0; idx < unvisitedCount; idx++) {
            int i = unvisited[idx];
            for (int j = 0; j < tourSize - 1; j++) {
                double cost = distanceMatrix[tour[j]][i] + distanceMatrix[i][tour[j + 1]] - distanceMatrix[tour[j]][tour[j + 1]];
                if (cost < localMinCosts[threadID]) {
                    localMinCosts[threadID] = cost;
                    localMinCostIndex[threadID] = i;
                    localInsertPositions[threadID] = j + 1;
                }
            }
        }
        #pragma omp critical
        {
           if (localMinCosts[threadID] < globalminCost) {
                    globalminCost = localMinCosts[threadID];
                    globalminCostIndex = localMinCostIndex[threadID];
                    globalinsertPosition = localInsertPositions[threadID];
                } 
        }
        }
        if (globalminCostIndex != -1) {
            for (int i = tourSize; i > globalinsertPosition; i--) {
                tour[i] = tour[i - 1];
            }
            tour[globalinsertPosition] = globalminCostIndex;

            for (int i = 0; i < unvisitedCount; i++) {
                if (unvisited[i] == globalminCostIndex) {
                    unvisited[i] = unvisited[unvisitedCount - 1];
                    unvisitedCount--;
                    break;
                }
            }

            tourSize++;
        }
    }
    
    free(unvisited);
    free(localMinCosts);
    free(localMinCostIndex);
    free(localInsertPositions);
    return tour;
}

 



int main(int argc, char *argv[]) {
    // Ensure correct usage
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    char *inputFile = argv[1];
    char *outputFile = argv[2];


    
    // Read coordinates
    int numOfCoords = readNumOfCoords(inputFile);
    double **coords = readCoords(inputFile, numOfCoords);
    
    // Generate distance matrix
    double **distanceMatrix = generateDistanceMatrix(coords, numOfCoords);

    // Apply the cheapest insertion algorithm
    int *tour = cheapestInsertion(distanceMatrix, numOfCoords);

    // Write the tour to the output file
    writeTourToFile(tour, numOfCoords + 1, outputFile);

    // Free memory
    for (int i = 0; i < numOfCoords; i++) {
        free(coords[i]);
        free(distanceMatrix[i]);

    
    }
    
    
    free(coords);
    free(distanceMatrix);
    free(tour);

    return 0;
}  