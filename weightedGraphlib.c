#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include"dllist.h"
#include "weightedGraphlib.h"
#include "fields.h"
#include <string.h>


#define INFINITIVE_VALUE 10000000
#define MAX_SIZE 100000


Graph createGraph(){
  Graph g = (Graph)malloc(sizeof(struct _Graph));
  if(g==NULL) {
	  printf("NULL, can not allocate mem, abort...\n");
	  return NULL;
  }
  g->vertices = make_jrb();
  g->edges = make_jrb();
  return g;
}

void addVertex(Graph graph, int id, char* name)
{
	JRB vertices  = graph->vertices;
	JRB edges  = graph->edges;
	JRB tree;
	
	//update vertex
	detailVertex iver = (detailVertex)malloc(sizeof(struct _detailVertex));
	if(iver==NULL) {
	  printf("NULL, can not allocate mem, abort...\n");
	  return;
	}
	iver->name = strdup(name);
	iver->indegreeTree = make_jrb();
	jrb_insert_int(vertices,id,new_jval_v(iver));
	
	//update edge
	tree = make_jrb();
	jrb_insert_int(edges,id,new_jval_v(tree));
}

char *getVertex(Graph graph, int id)
{
	JRB vnode=jrb_find_int(graph->vertices,id);
	
	if(vnode==NULL) return NULL;
	detailVertex iver = (detailVertex) vnode->val.v;
	
	return iver->name;
}

void addEdge(Graph graph, int v1, int v2, double wei)
{
	JRB enode = jrb_find_int(graph->edges,v1);
	JRB vnode = jrb_find_int(graph->vertices,v2);
	JRB tree;
	if((enode==NULL)||(vnode==NULL)) {
		printf("vertex not found\n");
		return;
    	}

	tree=(JRB)(enode->val).v;
	jrb_insert_int(tree,v2,new_jval_d(wei));
	
	detailVertex iver = (detailVertex) vnode->val.v;
	tree = iver->indegreeTree;
	jrb_insert_int(tree,v1,JNULL);
}

int hasEdge(Graph graph, int v1, int v2)
{
	JRB enode = jrb_find_int(graph->edges,v1);
	JRB tree;
	if(enode==NULL) {
		printf("vertex not found\n");
		return 0;
    	}
    	tree=(JRB)(enode->val).v;
	if(jrb_find_int(tree,v2)!=NULL)
		return 1;
	else return 0;
}

double getEdgeValue(Graph graph, int v1, int v2)
{
	JRB enode = jrb_find_int(graph->edges,v1);
	JRB tree;

	if(enode==NULL) {
		printf("vertex not found\n");
		return INFINITIVE_VALUE;
    	}
    
    	tree = jrb_find_int((JRB)(enode->val).v,v2);

	if(tree==NULL){
		return INFINITIVE_VALUE;
	}

    	return (tree->val).d;
}

int indegree(Graph graph, int v, int* output)
{
	JRB vnode = jrb_find_int(graph->vertices,v);
	JRB tree,node;
	int innum = 0;
	
	if(vnode==NULL) {
		printf("vertex not found\n");
		return 0;
    	}
    	detailVertex iver = (detailVertex) vnode->val.v;
	tree = iver->indegreeTree;
	//traverse tree
	jrb_traverse(node, tree){
		output[innum] = node->key.i;
		innum++;
		//printf("\nnumber innum:%d\n",output[innum-1]);
	}
	return innum;
}

int outdegree(Graph graph, int v, int* output)
{
	JRB enode = jrb_find_int(graph->edges,v);
	JRB tree,node;
	int outnum = 0;
	
	if(enode==NULL) {
		printf("vertex not found\n");
		return 0;
    	}

	tree =(JRB)(enode->val).v;
	//traverse tree
	jrb_traverse(node, tree){
		output[outnum] = node->key.i;
		outnum++;
		//printf("\nnumber innum:%d\n",output[outnum-1]);
	}
	return outnum;
}

void dropGraph(Graph graph)
{
	JRB node,innode;
	detailVertex detailnode;
	
	jrb_traverse(node,graph->edges){
		jrb_free_tree((JRB)jval_v(node->val));
	}
	jrb_free_tree(graph->edges);
	
	jrb_traverse(node,graph->vertices){
		detailnode = (detailVertex) node->val.v;
		free(detailnode->name);
		jrb_free_tree((JRB)detailnode->indegreeTree);
	}
	jrb_free_tree(graph->vertices);
	
	free(graph);
}

int extractMin(Dllist q, double setd[])
{
  	Dllist tmp, check;
  	int u;
  	double min = INFINITIVE_VALUE;
  	check = q->flink;
  	u = check->val.i;
  	dll_traverse(tmp, q)
  	{
   		if(setd[tmp->val.i] < min)
   		{
      			min = setd[tmp->val.i];
      			u = tmp->val.i;
      			check = tmp;
    		}
  	}
  	dll_delete_node(check);
  	return u;
}

void prim(Graph g, int s, int parent[])
{
  	double *cost = (double*)malloc(sizeof(double)*MAX_SIZE);
  	int *prev = (int*)malloc(sizeof(int)*MAX_SIZE);
  	int *output = (int*)malloc(sizeof(int)*MAX_SIZE);
  	int u, n, i;
  	Dllist queue = new_dllist();
  	JRB tmp;
  	jrb_traverse(tmp, g->vertices){
    		cost[tmp->key.i] = 10000;
    		prev[tmp->key.i] = -1;
    		dll_append(queue, new_jval_i(tmp->key.i));
  	}
  	cost[s] = 0;
  	while(!dll_empty(queue))
  	{
    		u = extractMin(queue, cost);
    		parent[u] = prev[u];
    	// printf("\npop %d : %d %d", u, cost[u], parent[u]);
   	 	n = outdegree(g, u, output);
    		for(i = 0; i < n; i++)
    		{
      			if(cost[output[i]] > getEdgeValue(g, u, output[i]))
      			{
        			cost[output[i]] = getEdgeValue(g, u, output[i]);
        			prev[output[i]] = u;
        // printf("\n%d -- > %d\n", u, output[i]);
      			}
    		}
  	}
  	free(prev);
  	free(cost);
  	free(output);
}

void printMST(Graph g, char* path, int size){
	int *parent = (int*)malloc(sizeof(int) * (size + 1));
  	FILE *p = fopen(path, "w");
  	if(p == NULL){
    		perror(path);
    		return;
  	}
	prim(g, 1, parent);
  	fprintf(p, "graph dothi");
  	fprintf(p, "\n{");
  	JRB tmp;
  	int *output = (int*)malloc(sizeof(int)*MAX_SIZE);
   	jrb_traverse(tmp, g->vertices){
     		if(-1 == parent[tmp->key.i]) continue;
    		fprintf(p, "\n%d -- %d [color = cyan]", tmp->key.i, parent[tmp->key.i] );
   	}
  	fprintf(p, "\n}");
  	free(output);

// JRB node;
// 	int* visited = (int*)calloc((size+1),sizeof(int));
// 	int* output = (int*)malloc(size*sizeof(int));
// 	jrb_traverse(node,g->vertices)
// 	{
// 		int n = outdegree(g,node->key.i,output);
// 		for(int i=0; i<n; i++)	
// 		{
// 			if(!visited[output[i]])
// 			{
// 				if(output[i] == parent[node->key.i] || node->key.i == parent[output[i]])
// 				fprintf(p,"\t%d -- %d [color = cyan];\n",node->key.i,output[i]);
// 				else
// 				fprintf(p,"\t%d -- %d;\n",node->key.i,output[i]);
// 			}
// 		}
// 		visited[node->key.i] = 1;
// 	}
//   	fprintf(p, "\n}");
// 	free(visited);
// 	free(output);

	free(parent);
	fclose(p);
}

void readFile(Graph g, char *filename, char *filedot,int*size) //doc file
{
  	IS is;
  	int numv = 0, nume = 0;
  	int i, n = 0;
  	is = new_inputstruct(filename);
  	if(is == NULL)
  	{
    		perror(filename);
    		return;
  	}
  	FILE *p = fopen(filedot, "w");
  	if(p == NULL)
  	{
    		printf("\nerror!\n");
    	return;
  	}
  	get_line(is);
  	char* type = is->fields[0];
  	if(strcmp(type, "[Weighted_graph]") == 0)
  	{
    		fprintf(p, "graph x");
    		fprintf(p, "\n{");
    		get_line(is);
    		numv = atoi(is->fields[0]);
    		nume = atoi(is->fields[1]);
    		for(i = 0 ; i < numv; i++)
    		{
      			get_line(is);
      			addVertex(g, atoi(is->fields[0]), is->fields[1]);
      			fprintf(p,"\n%d [style = filled, fillcolor = lightskyblue]", atoi(is->fields[0]));
    		}
    		for(i = 0; i < nume; i++)
    		{
      			get_line(is);
      			addEdge(g, atoi(is->fields[0]), atoi(is->fields[1]), atoi(is->fields[2]));
      			addEdge(g, atoi(is->fields[1]), atoi(is->fields[0]), atoi(is->fields[2]));
      			fprintf(p, "\n%d -- %d [label = %d ]", atoi(is->fields[0]), atoi(is->fields[1]), atoi(is->fields[2]));
    		}
    		fprintf(p, "\n}");
  	}
  	
  	if(strcmp(type, "[Directed_graph]") == 0)
  	{
    		fprintf(p, "digraph x");
    		fprintf(p, "\n{");
    		get_line(is);
    		numv = atoi(is->fields[0]);
    		nume = atoi(is->fields[1]);
    		for(i = 0 ; i < numv; i++)
    		{
      			get_line(is);
      			addVertex(g, atoi(is->fields[0]), is->fields[1]);
      			fprintf(p,"\n%d [style = filled, fillcolor = lightskyblue]", atoi(is->fields[0]));
    		}
    	for(i = 0; i < nume; i++)
    	{
      		get_line(is);
      		addEdge(g, atoi(is->fields[0]), atoi(is->fields[1]), 1);
      		fprintf(p, "\n%d -> %d", atoi(is->fields[0]), atoi(is->fields[1]));
    	}
    	fprintf(p, "\n}");
  	}

  	if(strcmp(type, "[Undirected_graph]") == 0)
  	{
    		fprintf(p, "graph x");
    		fprintf(p, "\n{");
    		get_line(is);
    		numv = atoi(is->fields[0]);
    		nume = atoi(is->fields[1]);
    		for(i = 0 ; i < numv; i++)
    		{
      			get_line(is);
      			addVertex(g, atoi(is->fields[0]), is->fields[1]);
      			fprintf(p,"\n%d [style = filled, fillcolor = lightskyblue]", atoi(is->fields[0]));
    		}
    		for(i = 0; i < nume; i++)
    		{
      			get_line(is);
      			addEdge(g, atoi(is->fields[0]), atoi(is->fields[1]), 1);
      			addEdge(g, atoi(is->fields[1]), atoi(is->fields[0]), 1);
      			fprintf(p, "\n%d -- %d", atoi(is->fields[0]), atoi(is->fields[1]));
    		}		
    		fprintf(p, "\n}");
  	}
  	fclose(p);
  	*size = numv;
  	jettison_inputstruct(is);
}

void relax(Graph g, int u, int v, double setd[], int setparent[])
{
  	if(setd[v] > setd[u] + getEdgeValue(g, u, v))
  	{
    		setd[v] = setd[u] + getEdgeValue(g, u, v);
    		setparent[v] = u;
  	}
}

double dijkstra(Graph g, int s, int t, int path[], int *length)
{
  	int u;
  	JRB tmp;
  	int *output = (int*)malloc(sizeof(int)*MAX_SIZE);
  	double *setd= (double*)malloc(sizeof(double)*MAX_SIZE);
  	int *setparent = (int*)malloc(sizeof(int)*MAX_SIZE);
  	Dllist pq = new_dllist();
  	jrb_traverse(tmp, g->vertices)
  	{
    		setd[tmp->key.i] = INFINITIVE_VALUE;
    		setparent[tmp->key.i] = -1;
    		dll_append(pq, new_jval_i(tmp->key.i));
  	}
  	setd[s] = 0;
  	while(!dll_empty(pq))
  	{
    		u = extractMin(pq, setd);
    		int n = outdegree(g, u, output);
    		int i;
    		for(i = 0; i < n; i++)
    		{
      			relax(g,u,output[i],setd,setparent);
    		}
  	}
  	if(setd[t] == INFINITIVE_VALUE)
  	{
    		free_dllist(pq);
    		return INFINITIVE_VALUE;
  	}
  	int count = 0;
  	int parent = t;
  	double weight = setd[t];
  // In duong di ngan nhat vaof path[]
  	while(parent != s)
  	{
    		path[count] = parent;
    		parent = setparent[parent];
    		count++;
  	}
  	path[count] = s;
  	*length = count;
  	free_dllist(pq);
  	free(setd);
  	free(output);
  	free(setparent);
  	return weight;
}

void BFS(Graph graph, int start, int stop,int*parent, int size)
{
	int *visited = (int*)malloc((size+1)*sizeof(int));
	int *output = (int*)malloc(size*sizeof(int));
   	for(int i=0; i<=size; i++)
	{
		visited[i] = 0;
	}
	int n, i, u, v;
  	Dllist node, queue;
   	queue = new_dllist();
   	dll_append(queue, new_jval_i(start));
   	while ( !dll_empty(queue) )
   	{
      		node = dll_first(queue);
      		u = jval_i(node->val);
      		dll_delete_node(node);
      		if (visited[u] != 1) 
      		{
          		visited[u] = 1;
				// printf("%d\t", u);
      	    		if ( u == stop)
        		{
				break;           
	  		}
          		n = outdegree(graph, u, output);
          		for (i=0; i<n; i++)
          		{
              			v = output[i];      							
				if (!visited[v])
              			{	
              				if(parent[v]==0) parent[v] = u;		
                 			dll_append(queue, new_jval_i(v));
	      			}
          		}
      		}
   	}      
	free(output);
	free(visited);	
}

void exch(int a[], int l, int r)
{
	int tmp = a[l];
	a[l] = a[r];
	a[r] = tmp; 
}

int Partition(int a[], int left , int right)
{
	int i = left; int j = right + 1; int pivot = a[left];
	while(i < j) 
	{	
		i++;
	 	while(i <= right && a[i] <= pivot){i++;}
    		j--;
		while(j >= left && a[j] > pivot){j--;}
		exch(a,i,j);
    	}
    	exch(a,i,j); exch(a,j,left);
    	return j;
}

void QuickSort(int a[], int left, int right)
{
	if(left < right)
	{
		int pivot = Partition(a,left,right);
		QuickSort(a,left,pivot - 1);
		QuickSort(a,pivot + 1,right);
	}
}

void Coloring(Graph g, int size, int Color[])
{
	int *output = (int*)malloc(sizeof(int)*size);
	for(int i=1; i<=size; i++)
	{
		Color[i] = 1;
		int n = outdegree(g,i,output);
		int *coutput = (int*)malloc(n*sizeof(int));
		for(int j=0; j<n; j++)
		{
			coutput[j] = Color[output[j]];
		}
		QuickSort(coutput,0,n-1);
		for(int j=0; j<n; j++)
		{
			if(Color[i] == coutput[j])
			{
				Color[i]++;
			}
		}
		free(coutput);
	}
	free(output);
}

void printColor(Graph g,char *filedot,int size)
{
	int *Color = (int*)malloc(sizeof(int)*(size+1));
	memset(Color, 0, (size+1)*sizeof(int));
	Coloring(g, size, Color);
	FILE* f2 = fopen(filedot,"w");
	fprintf(f2,"graph dothi\n{\n");
	for(int i=1; i<=size; i++)
	{
		if(Color[i]==1)
		{
			fprintf(f2,"%d [fillcolor=green, style=filled];\n",i);
		}
		if(Color[i]==2)
		{
			fprintf(f2,"%d [fillcolor=red, style=filled];\n",i);
		}
		if(Color[i]==3)
		{
			fprintf(f2,"%d [fillcolor=yellow, style=filled];\n",i);
		}
		if(Color[i]==4)
		{
			fprintf(f2,"%d [fillcolor=orange, style=filled];\n",i);
		}
		if(Color[i]==5)
		{
			fprintf(f2,"%d [fillcolor=pink, style=filled];\n",i);
		}
		if(Color[i]==6)
		{
			fprintf(f2,"%d [fillcolor=black, style=filled];\n",i);
		}
		if(Color[i]==7)
		{
			fprintf(f2,"%d [fillcolor=brown, style=filled];\n",i);
		}
		if(Color[i]==8)
		{
			fprintf(f2,"%d [fillcolor=gray, style=filled];\n",i);
		}
		if(Color[i]==9)
		{
			fprintf(f2,"%d [fillcolor=white, style=filled];\n",i);
		}
		if(Color[i]==10)
		{
			fprintf(f2,"%d [fillcolor=blue, style=filled];\n",i);
		}
	}
	JRB node;
	int* visited = (int*)calloc((size+1),sizeof(int));
	int* output = (int*)malloc(size*sizeof(int));
	jrb_traverse(node,g->vertices)
	{
		int n = outdegree(g,node->key.i,output);
		for(int i=0; i<n; i++)	
		{
			if(!visited[output[i]])
			{
				fprintf(f2,"\t%d -- %d;\n",node->key.i,output[i]);
			}
		}
		visited[node->key.i] = 1;
	}
	free(visited);
	free(output);
	fprintf(f2,"}\n");
	fclose(f2);
}

void explore(Graph g, int v, int* visited, int *clock)
{
	visited[v] = 1;
	// printf("%d\t", v);
	int *output = (int*)malloc(sizeof(int)*MAX_SIZE);
	int n;
	pre[v] = *clock; (*clock)++;
	n = outdegree(g, v, output);
	for(int i=0; i<n; i++)
	{
		if(!visited[output[i]])
		{
			explore(g,output[i],visited, clock);
		}
	}
	post[v] = *clock; (*clock)++;
}

void DFS(Graph g,int size)
{	
	int clock =1;
	int *visited = (int*)calloc(size+1,sizeof(int));
	JRB node;
	jrb_traverse(node,g->vertices)
	{
		if(!visited[node->key.i])
		{
			explore(g,node->key.i,visited, &clock);
		}
	}	
	free(visited);
}

Graph reverse(Graph g, int size)
{
	Graph h = createGraph();
	int *output = (int*)malloc(size*sizeof(int));
	JRB node;
	int n;
	jrb_traverse(node,g->vertices)
	{
		addVertex(h,node->key.i,getVertex(g,node->key.i));
	}
	jrb_traverse(node,g->vertices)
	{
		n = outdegree(g,node->key.i,output);
		for(int i=0; i<n; i++)
		{
			addEdge(h,output[i],node->key.i,1);
		}
	}
	free(output);
	return h;
}

int max(int* post,int size)
{
	int max = 1;
	for(int i=1; i<=size; i++)
	{
		if(post[i]>post[max]) max = i;
	}
	
	return max;
}

void explore1(Graph g, int v, int* visited, int* count)
{
	visited[v] = 1;
	int *output = (int*)malloc(sizeof(int)*MAX_SIZE);
	post[v] = 0;
	*count = *count + 1;
	int n;
	n = outdegree(g, v, output);
	for(int i=0; i<n; i++)
	{
		if(!visited[output[i]])
		{
			explore1(g,output[i],visited,count);
		}
	}
}

int soTPLTManh(Graph g, int size)
{
	Graph h = reverse(g,size);
	DFS(h,size);
	int v;
	int *visited = (int*)calloc(size+1,sizeof(int));
	int scc = 0; int count = 0;
	while(count<size)
	{
		v = max(post,size);
		scc++;
		explore1(g,v,visited,&count);
	}
	dropGraph(h);
	free(visited);	
	return scc;
}

int soTPLT(Graph g,int size)
{
	int clock = 1;
	int cc = 0;
	int *visited = (int*)calloc(size+1,sizeof(int));
	JRB node;
	jrb_traverse(node,g->vertices)
	{
		if(!visited[node->key.i])
		{
			cc = cc+1;
			explore(g,node->key.i,visited, &clock);
		}
	}	
	free(visited);
	return cc;

}

int IsDAG(Graph g,int size)
{
	DFS(g,size);
	JRB node;
	int n;
	int* output = (int*)malloc((size+1)*sizeof(int)); 
	jrb_traverse(node,g->vertices)
	{
		n = outdegree(g, node->key.i, output);	
		for(int i=0;i<n;i++){
			if(post[node->key.i]<post[output[i]])return 0;
		}
	}	
	free(output);
	return 1;
}

void  Topology(Graph g,int size)
{
	printf("\n");
	if(!IsDAG(g,size))
	{
		printf("Khong co trat tu Topo.\n");
		return;
	}
	Dllist Q, node;
	int* output = (int*)malloc((size+1)*sizeof(int)); 
	int* output1 = (int*)malloc((size+1)*sizeof(int)); 
	int n;
	printf("Thu tu Topo:\n");
	Q = new_dllist();
	for(int i = 1; i<=size; i++)
	{
		if(indegree(g,i,output)==0)
		{
			dll_append(Q, new_jval_i(i));
		}

	}
	for(int i=1; i<=size; i++)
	{
		output1[i]  = indegree(g,i,output);
	}
	while( !dll_empty(Q))
	{
		node = dll_first(Q);
		n = jval_i(node->val);
		dll_delete_node(node);
		// puts(getVertex(g,n));
		printf("%d\t", n);
		int m = outdegree(g,n,output);
		for(int i=0;  i<m; i++)
		{
			output1[output[i]]--;
			if(output1[output[i]]==0)
			{
				dll_append(Q, new_jval_i(output[i]));
			}
		}
	}
	free(output);
	free(output1);
	free_dllist(Q);
	printf("\n");
}
