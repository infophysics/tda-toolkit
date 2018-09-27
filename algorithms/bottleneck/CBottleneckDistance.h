/*
 * CBottleneckDistance.h
 *
 *  Created on: Feb 21, 2012
 *      Author: miro
 */

#ifndef CBOTTLENECKDISTANCE_H_
#define CBOTTLENECKDISTANCE_H_

#include <iostream>
#include <vector>


class Generator {
    public:
        double birth;
        double death;
};

class Edge {
    public:
        int vertex_1;
        int vertex_2;
        double weight;

        Edge(int v1, int v2, double w) {
            vertex_1 = v1;
            vertex_2 = v2;
            weight = w;
        }

        bool operator<(const Edge& other) const {
            return weight < other.weight;
        }
};

class CBottleneckDistance {
    private:

        bool NeglectShortGenerators;
        bool NeglectLastGenerators;

        double  NeglectSize;
        double  NeglectBornAfter;

        /**
         * Generators for two persistence diagrams which are
         * going to be compared
         */
        std::vector<Generator> Generators1;
        std::vector<Generator> Generators2;

        /**
         * Number of generators in Generator1 and Generator2
         * Total number of generators
         */
        unsigned int Max_Size;

        /**
         * Edges between all the nodes given by
         * generators of the persistence diagrams
         */
        std::vector<Edge> Edges;

        /**
         * Pairing between the vertices in the persistence
         * diagrams -1 means the vertex is not paired to any
         * other vertex it is pair to special vertex NIL
         * non negative number is the index of the vertex to
         * which the given vertex is paired
         */
        std::vector<int> Pair;

        /**
         * Connection matrix gives the edges between the
         * vertexes in the first diagram and the second
         */
        std::vector< std::vector<int> > Connections;

        /**
         * Layers used in Hopcroft-Karp algorithm the layer
         * information is required for a NIL (special) vertex
         * and all vertices in Generators1 Hence 0 slot is used
         * for NIL and all the vertices in Generators1 are
         * shifted by one. So to read layer of the vertex 0 we
         * access the Layers[1]
         */
        std::vector<int> Layers;

        /* Loads generators from the file */
        void LoadGeneratorsFromFile(const char* fileName,
                                    std::vector<Generator> &generators,
                                    double maxLevel);

        /* Computes distance between two generators */
        double  InfDistanceOfTwoGenerators(Generator gen1, Generator gen2);

        /* Computes distance from the diagonal for the given generator */
        double InfDistanceOfGeneratorFromDiagonal(Generator gen);

        /**
         * Produces all the edges between the nodes given by the
         * generators of the persistence diagrams stored in Generators1
         * and Generators2. The diagram are augmented to the same
         * length by projections to the diagonal.
         */
        void PrepareEdges();

        /* Depth first search used by Hopf-Karp algorithm */
        bool DFS(int v);

        /* Breath first search used by Hopf-Karp algorithm */
        bool BFS();

        /* Hopf-Karp algorithm to find maximal matching */
        void Hopcroft_Karp(unsigned int &matching);

    public:
        CBottleneckDistance();
        virtual ~CBottleneckDistance();

        /**
         *
         */
        double Distance(double maxLevel);

        /**
         *
         */
        double Distance(std::vector<Generator> generator1,
                        std::vector<Generator> generator2,
                        double maxLevel);
        /**
         *
         */
        double Distance(const char* diagram_1,
                        const char* diagram_2,
                        double maxLevel);

        /* Resets to default setings. It all generators are taken in account */
        void ResetSettings();

        /* Neglect the generatres woth live spam shorter than size */
        void SkipShortGenerators(double size);

        /* Neglect generators born aftre bornAfter */
        void SkipLastGenerators(double bornAfter);
};

#endif /* CBOTTLNECKDISTANCE_H_ */
