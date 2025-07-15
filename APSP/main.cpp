#include <iostream>
#include <boost/algorithm/string.hpp>  // Boost string utilities

#include "Dijkstra.h"
#include "algoritms/HiddenPaths.h"
#include "algoritms/UniformPaths.h"
#include "algoritms/Tree.h"
#include "algoritms/SCC.h"
#include "algoritms/Toroslu.h"
#include "algoritms/DijkstraAPSP.h"
#include "graphGenerators.h"
#include "utils.h"
#include "run_algorithm.h"
#include "GraphTools/ba_graph.h"

using namespace std;

int main() {
    srand(time(NULL));
    // UniformPathsAlgorithm, HiddenPathsAlgorithm, floydWarshall, DijkstraAPSPc, Toro
    // "Uniform", "Hidden", "FW", "Dijkstra", "Toroslu"
    AlgorithmFP algs[] = {SCCalgo, FastHourglassHalf, floydWarshall, DijkstraAPSPc, Toro, UniformPathsAlgorithm, HiddenPathsAlgorithm, NULL}; // vedno dodaj null na koncu!
    char* algNames[] = {"SCC", "Tree", "FW", "Dijkstra", "Toroslu", "Uniform", "Hidden", "NULL"};

    RunTurkTest(algs, algNames);      // ErdosRenyi
    //RunTurk2Test(algs, algNames);     // TurkTest with extra data (sccV, sccE, wccV, wccE)
    //RunN_1Test(algs, algNames);       // CompGraph
    //RunlgNTest(algs, algNames);       // LogGraph
    //RunSparseTest(algs, algNames);    // Cycle + ErdosRenyi
    //RunSCCsTest(algs, algNames);
    //RunCompleteTest(algs, algNames, 0);
    //RunGnutellaTest(algs, algNames);
    //RunTxtTest(algs, algNames);
    //Run_BA_Test(algs, algNames);

    /*int n = 8;
    float** graph = generate_BA_graph(n, log(n));
    print2DGraph(graph, n);
    freeGraph(graph, n);*/

    return 0;
}
