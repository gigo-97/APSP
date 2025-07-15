//
// Created by Dani Zugan on 28. 5. 25.
//

#ifndef BA_GRAPH_H
#define BA_GRAPH_H

#include "../utils.h"
#include "../run_algorithm.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

//float** generate_BA_graph(int n, int m);
void Run_BA_Test(AlgorithmFP* algs, char** algNames);

#ifdef __cplusplus
}
#endif

#endif //BA_GRAPH_H
