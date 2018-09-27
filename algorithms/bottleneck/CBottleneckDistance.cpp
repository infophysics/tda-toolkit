/**
 * CBottleneckDistance.cpp
 *
 *  Created on: Feb 21, 2012
 *      Author: miro
 *      Modified by: Kelvin Abrokwa (kelvinabrokwa@gmail.com)
 */

#include "CBottleneckDistance.h"

#include <fstream>
#include <math.h>
#include <queue>
#include <algorithm>


/**
 * Debugging methods
 */
void printGenerators(std::vector<Generator> generators)
{
    for (auto gen : generators)
        std::cout << gen.birth << " - " << gen.death << std::endl;
}

CBottleneckDistance::CBottleneckDistance()
{
    NeglectShortGenerators = false;
    NeglectLastGenerators = false;
}

/**
 * Destructor
 */
CBottleneckDistance::~CBottleneckDistance() {}

/**
 *
 */
void CBottleneckDistance::LoadGeneratorsFromFile(
        const char* fileName, std::vector<Generator> &generators, double maxLevel)
{

    Generator gen;

    // erase generators
    generators.clear();

    // Open the file
    std::ifstream infile(fileName,  std::ifstream::in);

    //Test if the file was sucssefuly opened
    if (!infile) {
        std::cout
            << " Unable to open input file "
            << fileName
            << ". I'll work with an empty persistence diagram!"
            << std::endl;
        return;
    }

    // Read the generators from the file
    while (!infile.eof()) {
        // reads one line of the file
        infile >> gen.birth;

        // Just in case that files ends in a strange way
        if (infile.eof())
            break;
        infile >> gen.death;

        // Encorporate Vidits standard that the generator which lives till the end dies at -1
        if (gen.birth == -1)
            gen.birth = maxLevel;
        if (gen.death == -1)
            gen.death = maxLevel;

        // Test if generator should be included (length of life spam plus birth)
        if (NeglectShortGenerators && (gen.death - gen.birth <= NeglectSize))
            continue;
        if (NeglectLastGenerators && (gen.birth >=  NeglectBornAfter))
            continue;

        generators.push_back(gen);
    }

    // Close the file
    infile.close();

    return;
}

/**
 *
 */
double CBottleneckDistance::InfDistanceOfTwoGenerators(Generator gen1, Generator gen2)
{
    double birth_diff = fabs(gen1.birth - gen2.birth);
    double death_diff = fabs(gen1.death - gen2.death);

    if (birth_diff > death_diff)
        return 0; // birth_diff;

    return 0; // death_diff;
}


/**
 *
 */
double CBottleneckDistance::InfDistanceOfGeneratorFromDiagonal(Generator gen)
{
    return fabs(gen.death - gen.birth) / 2; // fabs( gen.death - gen.birth) / 2.0;
}

/**
 * Prepare Edges
 */
void CBottleneckDistance::PrepareEdges()
{
    // Set the number of generators
    Max_Size = Generators1.size() + Generators2.size();

    // Clear edges
    Edges.clear();

    //	std::cout << "Max_size = " << Max_Size                      << "\n";
    //	std::cout << "n gen 1  = " << Generators1.size()            << "\n";
    //	std::cout << "n diag 1 = " << Max_Size - Generators1.size() << "\n";
    //	std::cout << "n gen 2  = " << Generators2.size()            << "\n";
    //	std::cout << "n diag 2 = " << Max_Size -Generators2.size()  << "\n";

    // Connect all diagonal points to each other
    for (unsigned int i = Generators1.size(); i < Max_Size; i++)
        for (unsigned int j = Max_Size + Generators2.size(); j < 2 * Max_Size; j++)
            Edges.push_back(Edge(i, j, 0));

    //  int number_of_pure_diag_edges = Edges.size();
    //  std::cout << "Number of pure diag edges = " << number_of_pure_diag_edges <<"\n";

    // Edges between real points
    unsigned int i = 0;
    for (std::vector<Generator>::const_iterator cur1 = Generators1.begin(); cur1 != Generators1.end(); cur1++) {
        unsigned int j = Max_Size;
        for (std::vector<Generator>::const_iterator cur2 = Generators2.begin(); cur2 != Generators2.end(); cur2++) {
            Edges.push_back(Edge(i, j, InfDistanceOfTwoGenerators(*cur1, *cur2)));
            j++;
        }
        i++;
    }

    //	 int number_of_real_points_edges = Edges.size() - number_of_pure_diag_edges;
    //	 int sum2 = Edges.size();
    //	 std::cout << "Number of real edges = " << number_of_real_points_edges <<"\n";

    // Edges between real points and their corresponding diagonal points
    i = 0;
    for (std::vector<Generator>::const_iterator cur1 = Generators1.begin(); cur1 != Generators1.end(); ++cur1, ++i)
        Edges.push_back( Edge(i, Max_Size + Generators2.size() + i, InfDistanceOfGeneratorFromDiagonal(*cur1)));

    //	int nn = Edges.size() - sum2;
    //	std::cout << "nn = " << nn <<"\n";

    i = Max_Size;
    for (std::vector<Generator>::const_iterator cur2 = Generators2.begin(); cur2 != Generators2.end(); ++cur2, ++i)
        Edges.push_back( Edge( Generators1.size() + (i - Max_Size), i, InfDistanceOfGeneratorFromDiagonal(*cur2)));

    std::sort(Edges.begin(), Edges.end());

    //	std::cout << "N edges " << Edges.size() << "\n";
    //	for( i = 0; i < Edges.size(); ++i)
    //		std::cout
    //		    << i
    //		    << " -------- "
    //		    << Edges[i].vertex_1
    //		    << " "
    //		    << Edges[i].vertex_2
    //		    << " "
    //		    << Edges[i].weight
    //		    << std::endl;
}


/*
 * Depth first search used by Hopf-Karp algorithm
 */
bool CBottleneckDistance::DFS(int v)
{
    if (v < 0)
        return true;

    // for every adjacent vertex u of v
    for (unsigned int i = 0; i < Connections[v].size(); ++i) {
        int u = Connections[v][i];
        if (Layers[Pair[u] + 1] == Layers[v + 1] + 1) {
            if (DFS( Pair[u])) {
                Pair[u] = v;
                Pair[v] = u;
                return true;
            }
        }
    }

    Layers[v + 1] = -1;

    return false;
}


/**
 * Breath first search used by Hopcroft-Karp algorithm
 */
bool CBottleneckDistance::BFS()
{
    std::queue<int> vertex_queue;

    // For every vertex v given by Generators1
    for (unsigned int v = 0; v < Max_Size; v++) {
        // If its not paired to vertex in Generators2
        if (Pair[v] < 0) {
            // Set its layer to 0 and put it in the queue
            Layers[v + 1] = 0;
            vertex_queue.push(v);
        }
        else {
            // Otherwise mark it as matched (set Layer to NILL)
            Layers[v + 1] = -1;
        }
    }

    // Set layer for NIL
    Layers[0] = -1;

    // Search the vertices in the queue
    while (!vertex_queue.empty()) {
        int v = vertex_queue.front();
        vertex_queue.pop();
        if (Layers[ v + 1 ] > Layers[ 0 ]) {
            for(unsigned int i = 0; i < Connections[ v ].size(); ++i) {
                int u = Connections[ v ][ i ];
                // Check if the vertex has an edge to the match vertex
                if (Layers[ Pair[ u ] + 1 ] < 0) {
                    // Set the layer of the vertex (it can be NILL) which is matched to the matched vertex u
                    Layers[ Pair[ u ] + 1 ] = Layers[ v + 1 ] + 1;
                    // If the pairing vertex is not NIL add it into the queue
                    if(Pair[ u ] != -1)
                        vertex_queue.push(Pair[u]);
                }
            }
        }
    }

    return Layers[0] != -1;
}


/**
 * Hopf-Karp algorithm to find maximal matching
 */
void CBottleneckDistance::Hopcroft_Karp(unsigned int &matching)
{
    while (BFS() == true)
        for (unsigned int vertex = 0; vertex < Max_Size; ++vertex) {
            if (Pair[ vertex ] == -1) {
                if (DFS(vertex)) {
                    ++matching;
                }
            }
        }
}


/**
 * A ::Distance wrapper that takes a vector of Generators
 */
double CBottleneckDistance::Distance(
        std::vector<Generator> generators1, std::vector<Generator> generators2, double maxLevel)
{
    Generators1 = generators1;
    Generators2 = generators2;
    return Distance(maxLevel);
}


/**
 * A ::Distance wrapper that loads persistence diagrams from a file
 */
double CBottleneckDistance::Distance(const char* diagram_1, const char* diagram_2, double maxLevel)
{
    // Load the diagrams
    LoadGeneratorsFromFile(diagram_1, Generators1, maxLevel);
    LoadGeneratorsFromFile(diagram_2, Generators2, maxLevel);
    return Distance(maxLevel);
}

/**
 * Calculate bottleneck distance of Generators1 and Generators2
 */
double CBottleneckDistance::Distance(double maxLevel)
{
    // If both diagrams are empty the distance is 0
    if (Generators1.size() == 0 &&  Generators2.size() == 0)
        return 0;

    PrepareEdges();

    // Clear the pairing
    Pair.clear();
    Pair.assign(2 * Max_Size, -1);

    // Clearing Layers
    Layers.clear();
    Layers.resize(Max_Size + 1);

    // No vertices are matched
    unsigned int matching = 0;

    // Clear the connection matrix and set it to the right size
    Connections.clear();
    Connections.resize(Max_Size);

    // The maximal weight of the edges which are used for the matching
    double current_weight = 0;

    // First non added edge is an iterator pointing to the first edge
    // in Edges which was added to the Connections */
    unsigned int first_not_added_edge = 0;

    // Repeat till all the vertices are matched
    while(matching < Max_Size) {

        // Add the edges with the current weight (distance) to the connection matrix
        while (Edges[first_not_added_edge].weight == current_weight && first_not_added_edge < Edges.size()) {
            // Add the edge to Connections
            Connections[Edges[first_not_added_edge].vertex_1].push_back(Edges[first_not_added_edge].vertex_2);
            ++first_not_added_edge;
        }

        // Do matching
        Hopcroft_Karp(matching);
        if (matching == Max_Size)
            return current_weight;

        // Check if we did not run out of edges. This should never happen.
        if (first_not_added_edge == Edges.size()) {
            std::cout << "Serious problem - Not enough eddges to find the matching! \n";
            return -1;
        }

        // Increase the value of the current weight
        current_weight = Edges[ first_not_added_edge ].weight;
    }

    return -1;
}


/**
 *
 */
void CBottleneckDistance::ResetSettings()
{
    NeglectShortGenerators = false;
    NeglectLastGenerators = false;
}


/**
 *
 */
void CBottleneckDistance::SkipShortGenerators(double size)
{
    NeglectShortGenerators = true;
    NeglectSize = size;
}


/**
 *
 */
void CBottleneckDistance::SkipLastGenerators(double bornAfter)
{
    NeglectLastGenerators = true;
    NeglectBornAfter = bornAfter;
}

