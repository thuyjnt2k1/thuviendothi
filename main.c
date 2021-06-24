#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include<string.h>
#include<time.h>
#include "weightedGraphlib.h"

void makeWeightedGraph(char *pathname,int num){
	FILE *p = fopen(pathname, "w");
	if(p == NULL){
		perror(pathname);
		return;
	}
	int u;
	int count = 1;
	fprintf(p, "[Weighted_graph]\n");
	srand((int) time(0));
	int n;
	do{
		n = rand()%(num*2);
	}while(n < num);
	fprintf(p, "%d %d", num, n);
	int i = 0;
	for(i = 1; i <= num; i++){
		fprintf(p, "\n%d %d", i, i);
	}
	for(i = 0; i < n; i++){
		fprintf(p, "\n%d %d %d",rand()%num+1, rand()%num+1, rand()%100 );
	}


	fclose(p);
}

void makeDirectedGraph(char *pathname,int num){
	FILE *p = fopen(pathname, "w");
	if(p == NULL){
		perror(pathname);
		return;
	}
	int u;
	int count = 1;
	fprintf(p, "[Directed_graph]\n");
	srand((int) time(0));
	int n;
	do{
		n = rand()%(num*2);
	}while(n < num);
	fprintf(p, "%d %d", num, n);
	int i = 0;
	for(i = 1; i <= num; i++){
		fprintf(p, "\n%d %d", i, i);
	}
	for(i = 0; i < n; i++){
		fprintf(p, "\n%d %d",rand()%num+1, rand()%num+1);
	}
	fclose(p);
}

void makeUnDirectedGraph(char *pathname,int num){
	FILE *p = fopen(pathname, "w");
	if(p == NULL){
		perror(pathname);
		return;
	}
	int u;
	int count = 1;
	fprintf(p, "[Undirected_graph]\n");
	srand((int) time(0));
	int n;
	do{
		n = rand()%(num*2);
	}while(n < num);
	fprintf(p, "%d %d", num, n);
	int i = 0;
	for(i = 1; i <= num; i++){
		fprintf(p, "\n%d %d", i, i);
	}
	for(i = 0; i < n; i++){
		fprintf(p, "\n%d %d",rand()%num+1, rand()%num+1);
	}
	fclose(p);
}

int main()
{
	int i;
	pre = (int*)malloc(sizeof(int) * MAX_SIZE);
	post = (int*)malloc(sizeof(int) * MAX_SIZE);
	int *output = (int*)malloc(sizeof(int) * MAX_SIZE);
	memset(output, 0, MAX_SIZE*sizeof(int));
  	Graph g = createGraph();
  	int size;
	makeUnDirectedGraph("text3.txt", 50);
	makeDirectedGraph("text1.txt", 50);
	makeWeightedGraph("text2.txt", 50);
	readFile(g, "text3.txt", "graph.dot", &size);

	BFS(g, 1, -1, output, size);
	int i = 10;
	while(i != 1)
	{
	 	printf("%d\t", i);
	 	i = output[i];
	}
	printf("%d\t", i);

//	printColor(g,"color.dot",size);

//	printf("%d", soTPLTManh(g, size));

//	printf("\n%d",soTPLT(g,size));
	
	// int n;
	// printf("tong quang duong: %.0lf\nquang duong: ", dijkstra(g, 1, 25, output, &n));
	// for(i = n; i >= 0; i--){
	// 	printf("%d ", output[i]); 
	// }

	// Topology(g, size);

	// printMST(g, "prim.dot", size);

	// DFS(g,size);
	
	// for(j=1;j<=size;j++){
	// 	printf("%d %d ", pre[j], post[j]);
	// 	printf("\n");
	// }

	// if(IsDAG(g,size)){
	// 	printf("\n DAG !");
	// }
	// else printf("\n Not DAG");

	free(pre);
	free(post);
	free(output);
	dropGraph(g);
	printf("\n");
	return 0;
}
