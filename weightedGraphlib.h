#ifndef weightedGraphlib_h
#define weightedGraphlib_h
#endif

#define INFINITIVE_VALUE 10000000
#define MAX_SIZE 100000


#include "jrb.h"

int*pre;
int*post;

typedef struct _detailVertex{
	char *name;
	JRB indegreeTree;
} *detailVertex;

typedef struct _Graph{
	JRB edges;
	JRB vertices;
} *Graph;

Graph createGraph();
void addVertex(Graph graph, int id, char* name);
char *getVertex(Graph graph, int id);
void addEdge(Graph graph, int v1, int v2, double weight);
double getEdgeValue(Graph graph, int v1, int v2);
int hasEdge(Graph graph, int v1, int v2);
int indegree(Graph graph, int v, int* output);
int outdegree(Graph graph, int v, int* output);
void prim(Graph g, int s, int parent[]);
void printMST(Graph g, char* path, int size);
void readFile(Graph g, char *filename, char *filedot,int*size);
double dijkstra(Graph g, int s, int t, int path[], int *length);
void BFS(Graph graph, int start, int stop,int*parent, int size);
void printColor(Graph g,char *filedot,int size);
int soTPLTManh(Graph g, int size);
int soTPLT(Graph g,int size);
int IsDAG(Graph g,int size);
void DFS(Graph g,int size);
void  Topology(Graph g,int size);
void dropGraph(Graph graph);


