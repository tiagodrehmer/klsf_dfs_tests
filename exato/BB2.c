#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>


// gcc -o p BB2.c -lm


// #define MAX_ITERATION 1000 //1017
// #define MAX_IN_SOLUTION 10 //1017
// #define FRAC_FIX 0.8
// #define FRAC_CORES_SUB_SAI 0.2
// #define FRAC_SUB_FICA 0.8
// #define MAX_POND 1
//numero 2

int MAX_ITERATION; //1017
int MAX_IN_SOLUTION; //1017
int MAX_POND;
float FRAC_FIX;
float FRAC_CORES_SUB_SAI;
float FRAC_SUB_FICA;

typedef struct Edges{
	int v1, v2, l;
}Edges;

typedef struct Labels{
	int my_value;
	int qtd_edges;
	int *qtd_edge_per_k;
	Edges *es;
}Labels;




typedef struct Vertices{
	int qtd_adj;
	Edges *es;
}Vertices;

typedef struct Problem{
	int V;
	int L;
	int E;
	int k;
	Labels *root_labels;
	Labels *edges_per_label;
	Vertices *edges_per_vertice;
} Problem;

typedef struct COMPS{
	int qtd_vs;
	int *index_vs;
}COMPS;







typedef struct BB{
	COMPS *comps;
	int *visitados;
	Problem *subs;
	int **m_adj;
	int *qtd_edges_added;
	int *Vs;
	int *solution;

}BB;

typedef struct VNS{
	int *best_know;
	int *solution_now;
	int *ponderamento;
	int *labels_sub;
	int *labels_fixed;
	int *bool_labels;
	int *bool_solution;
	int iteracao;
	int qtd_labels_sub;
	int qtd_fixed;
	int total_pond;
	int max_in_solution;
	int best;
	int solution_value;
	long time;
	BB bb;
	int it_in_solution;
}VNS;


void 
input(){
	int x;
	printf("\n");
	scanf("%d", &x);
}

void 
printa_infos(VNS vns, int V, int E, int L, int k){
	//printf("\nqtd_labels_sub: %d", vns.qtd_labels_sub);
	printf("\nV: %d E: %d L:%d ", V, E, L);
	// printf("\nlabels_sub: ");
	// for(int i = 0; i < vns.qtd_labels_sub; i++)
	// 	printf("%d ", vns.labels_sub[i]);
	// printf("\nsolution_now: ");
	// for(int i = 0; i < k; i++)
	// 	printf("%d ", vns.solution_now[i]);
	// printf("\nbest_know: ");
	// for(int i = 0; i < k; i++)
	// 	printf("%d ", vns.best_know[i]);
	//printf("\nbool_labels: ");
	//for(int i = 0; i < L; i++)
	//	printf("%d ", vns.bool_labels[i]);
	//printf("\nponderamento: ");
	//for(int i = 0; i < L; i++)
	//	printf("%d ", vns.ponderamento[i]);

}


void 
printa_problem(Problem P){
	printf("\nV: %d, E: %d L: %d, k: %d", P.V, P.E, P.L, P.k);
	Edges e;
	for (int l = 0; l < P.L; l++){
		printf("\n%d: ", l);
		for(int i = 0; i < P.edges_per_label[l].qtd_edges; i++){
			e = P.edges_per_label[l].es[i];
			printf("(%d, %d)", e.v1, e.v2);
		}
	}
	printf("\n");
}



void 
add_edge(Problem *G, int v1, int v2, int l){			
	Edges e;
	e.v1 = v1;
	e.v2 = v2;
	e.l = l;
	// add edge per label
	G->edges_per_label[l].qtd_edges += 1;
	if (G->edges_per_label[l].qtd_edges == 1){
		G->edges_per_label[l].es = (Edges*) malloc(G->edges_per_label[l].qtd_edges * sizeof(Edges));
	}else{
		G->edges_per_label[l].es = (Edges*) realloc(G->edges_per_label[l].es, G->edges_per_label[l].qtd_edges * sizeof(Edges));
	}
	if(G->edges_per_label[l].es == NULL){
		printf("out of memory");
		exit(1);
	}
	G->edges_per_label[l].es[G->edges_per_label[l].qtd_edges - 1] = e;



	//add em v1
	G->edges_per_vertice[v1].qtd_adj += 1;

	if (G->edges_per_vertice[v1].qtd_adj == 1){
		G->edges_per_vertice[v1].es = (Edges*) malloc(G->edges_per_vertice[v1].qtd_adj * sizeof(Edges));
	}else{
		G->edges_per_vertice[v1].es = (Edges*) realloc(G->edges_per_vertice[v1].es, G->edges_per_vertice[v1].qtd_adj * sizeof(Edges));

	}

	if(G->edges_per_vertice[v1].es == NULL){
		printf("out of memory");
		exit(1);
	}
	
	G->edges_per_vertice[v1].es[G->edges_per_vertice[v1].qtd_adj - 1] = e;


	//add em v2
	G->edges_per_vertice[v2].qtd_adj += 1;

	if (G->edges_per_vertice[v2].qtd_adj == 1){
		G->edges_per_vertice[v2].es = (Edges*) malloc(G->edges_per_vertice[v2].qtd_adj * sizeof(Edges));
	}else{
		G->edges_per_vertice[v2].es = (Edges*) realloc(G->edges_per_vertice[v2].es, G->edges_per_vertice[v2].qtd_adj * sizeof(Edges));
	}

	if(G->edges_per_vertice[v2].es == NULL){
		printf("out of memory");
		exit(1);
	}
	
	G->edges_per_vertice[v2].es[G->edges_per_vertice[v2].qtd_adj - 1] = e;
}



void 
init_problem(Problem *G){
	G->edges_per_vertice = (Vertices*) malloc(G->V*sizeof(Vertices)); 
	if(G->edges_per_vertice == NULL){
		printf("Out of memory");
		exit(1);
	}
	for (int i = 0; i < G->V; i++){
		G->edges_per_vertice[i].es = NULL;
		G->edges_per_vertice[i].qtd_adj = 0;
	}
	int j;
	G->root_labels = NULL;
	G->edges_per_label = (Labels*) malloc(G->L * sizeof(Labels));
	for (int i = 0; i < G->L; i++){
		G->edges_per_label[i].qtd_edge_per_k = malloc((G->k + 1) * sizeof(int));
		for(j = 0; j < G->k + 1; j++){
			G->edges_per_label[i].qtd_edge_per_k[j] = 0;
		}
		G->edges_per_label[i].es = NULL;
		G->edges_per_label[i].qtd_edges = 0;
	}
}


int 
read_file(Problem *G, char* s_nome_arquivo){

	int v1, v2, l;




	FILE *f_arquivo_exec = fopen(s_nome_arquivo, "r");

	if (f_arquivo_exec == NULL){
		printf("Erro no arquivo");
		return(1);
	}
	


	fscanf(f_arquivo_exec, "%d %d %d %d", &G->V, &G->E, &G->L, &G->k);


	init_problem(G);
	int id = 0;
	while(fscanf(f_arquivo_exec, "%d %d %d", &v1, &v2, &l) != EOF){
		add_edge(G, v1-1, v2-1, l-1);
		id++;
	}

	int j;
	for (int i = 0; i < G->L; i++){
		G->edges_per_label[i].my_value = i;
		for(j = 0; j < G->k + 1; j++){
			G->edges_per_label[i].qtd_edge_per_k[j] = G->edges_per_label[i].qtd_edges;
		}
	}
	fclose(f_arquivo_exec);
	return(0);
}


Problem 
start_problem(){
	Problem P;
	P.V = 0;
	P.E = 0;
	P.L = 0;
	P.k = 0;
	P.edges_per_label = NULL;
	P.edges_per_vertice = NULL;
	return P;
}

void 
finish_vns(VNS *vns){
	if(vns->best_know != NULL)
		free(vns->best_know);
	if(vns->ponderamento != NULL)
		free(vns->ponderamento);
	if(vns->labels_fixed != NULL)
		free(vns->labels_fixed);
	if(vns->labels_sub != NULL)
		free(vns->labels_sub);
	if(vns->solution_now != NULL)
		free(vns->solution_now);
	if(vns->bool_labels != NULL)
		free(vns->bool_labels);
	if(vns->bool_solution != NULL)
		free(vns->bool_solution);
}

void 
finish_problem(Problem *G){
	if(G->edges_per_label != NULL){
		for (int l = 0; l < G->L; l ++){
			if(G->edges_per_label[l].es != NULL)
				free(G->edges_per_label[l].es);
			free(G->edges_per_label[l].qtd_edge_per_k);
		}
		free(G->edges_per_label);
	}
	if(G->edges_per_vertice != NULL){
		for (int v = 0; v < G->V; v ++){
			if(G->edges_per_vertice[v].es != NULL)
				free(G->edges_per_vertice[v].es);
		}
		free(G->edges_per_vertice);

	}
}



int 
DFS_avalia(int v, int componente, int *solution, int *visitados, Problem G){
	visitados[v] = componente;
	Edges e;
	for(int i = 0; i < G.edges_per_vertice[v].qtd_adj; i++){
		e = G.edges_per_vertice[v].es[i];
		if(solution[e.l] > 0.9){
			if(e.v1 == v){
				if(visitados[e.v2] == -1){
					DFS_avalia(e.v2, componente, solution, visitados, G);
				}
			}else{
				if(visitados[e.v1] == -1){
					DFS_avalia(e.v1, componente, solution, visitados, G);
				}
			}
		}
	}
}

int 
avalia_solution(int *solution, Problem G){
	float componente = 0;

	int visitados[G.V];

	for(int i = 0; i < G.V; i++)
		visitados[i] = -1;

	for(int v = 0; v < G.V; v++){
		if(visitados[v] == -1){
			componente = componente + 1;
			DFS_avalia(v, componente, solution, visitados, G);
		}
	}

	return(componente);
}



void to_file(VNS vns, int k, int L){
	FILE *saida = fopen("arquivo_aux.txt", "w");
	for(int i = 0; i < k; i++)
		fprintf(saida, "%d ", vns.solution_now[i]);
	fprintf(saida, "\n");
	for(int i = 0; i < k; i++)
		fprintf(saida, "%d ", vns.best_know[i]);
	fprintf(saida, "\n");
	for(int i = 0; i < L; i++)
		fprintf(saida, "%d ", vns.ponderamento[i]);
	fprintf(saida, "\n%d %d %ld %d %d %d", vns.best, vns.solution_value, vns.time, vns.iteracao, vns.it_in_solution, vns.total_pond);
	fclose(saida);
}




void merge(Labels arr[], int l, int m, int r, int k_atual)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
  
    /* create temp arrays */
    Labels L[n1], R[n2];
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i].qtd_edge_per_k[k_atual] >= R[j].qtd_edge_per_k[k_atual]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
  
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
  
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
  
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(Labels arr[], int l, int r, int k_atual)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(arr, l, m, k_atual);
        mergeSort(arr, m + 1, r, k_atual);
  
        merge(arr, l, m, r, k_atual);
    }
}

VNS 
start_alg(int V, int L, int k){
	VNS vns;
	vns.max_in_solution = MAX_IN_SOLUTION;
	vns.iteracao = 0;
	vns.ponderamento = (int*) malloc(L * sizeof(int));
	vns.best_know = (int*) malloc(k*sizeof(int));
	vns.bool_labels = (int*) malloc(L*sizeof(int));
	vns.qtd_fixed = (int) floor(k * FRAC_FIX);
	vns.qtd_labels_sub = k - vns.qtd_fixed;
	if(k == vns.qtd_fixed){
		vns.qtd_fixed --;
	}
	vns.labels_fixed = (int*) malloc(vns.qtd_fixed * sizeof(int));
	vns.labels_sub = (int*) malloc(L * sizeof(int));
	vns.solution_now = (int*) malloc(k*sizeof(int));
	vns.bool_solution = (int*) malloc(k*sizeof(int));
	vns.best = V;
	vns.solution_value = V;
	vns.time = time(NULL);
	vns.it_in_solution = 0;


	vns.total_pond = 0;
	for(int i = 0; i < L ; i++){
		vns.labels_sub[i] = -1;
	}

	for(int i = 0; i < vns.qtd_fixed; i++){
		vns.labels_fixed[i] = 0;
	}

	for(int i = 0; i < L; i++){
		vns.ponderamento[i] = MAX_POND;
		vns.bool_labels[i] = vns.iteracao;
	}
	vns.total_pond = k * MAX_POND;

	int cont_l = 0, r;

	while(cont_l < k){
		r = rand()%L;
		while(vns.bool_labels[r] != vns.iteracao){
			r = rand()%L;
		}
		vns.bool_labels[r] = vns.iteracao;
		vns.solution_now[cont_l] = r;
		vns.bool_solution[cont_l] = -1;
		vns.best_know[cont_l] = r;
		cont_l += 1;
	}
		
	
	return(vns);
}





COMPS *init_comps(int V, int* visitados){
	COMPS* comps = (COMPS*) malloc(V * sizeof(COMPS));
	for(int v = 0; v < V; v++){
		visitados[v] = v;
		comps[v].index_vs = (int*) malloc(V * sizeof(int));
		comps[v].qtd_vs = 1;
		for(int i = 0; i < V; i++){
			comps[v].index_vs[i] = -1;
		}
		comps[v].index_vs[0] = v;			
	}
	return(comps);
}

void finish_comps(COMPS *comps, int V){
	if(comps != NULL){
		for(int v = 0; v < V; v++){
			if(comps[v].index_vs != NULL)
				free(comps[v].index_vs);
			comps[v].qtd_vs = 0;
		}
		free(comps);
	}
}

void printa_comps(COMPS *comps, int V){
	for(int i = 0; i < V; i++){
		printf("\n %d: ", i);
		for(int j = 0; j < V; j++){
			printf("%d ", comps[i].index_vs[j]);
		}
	}
}


void fix_labels(int* visitados, COMPS* comps, Labels label, int *V){
	int v1, v2, aux;
	for(int i = 0; i < label.qtd_edges; i++){
		v1 = label.es[i].v1;
		v2 = label.es[i].v2;
		if(visitados[v1] != visitados[v2]){
			
			if(visitados[v1] < visitados[v2]){
				aux = visitados[v1];
				v2 = visitados[v2];
				v1 = aux;
			}else if(visitados[v2] < visitados[v1]){
				aux = visitados[v2];
				v2 = visitados[v1];
				v1 = aux;
			}

			while(comps[v2].qtd_vs){
				comps[v2].qtd_vs--;
				aux = comps[v2].index_vs[comps[v2].qtd_vs];
				comps[v1].index_vs[comps[v1].qtd_vs++] = aux;
				visitados[aux] = v1;
				comps[v2].index_vs[comps[v2].qtd_vs] = -1;
			}
			*V -= 1;
			while(comps[*V].qtd_vs){
				comps[*V].qtd_vs--;
				aux = comps[*V].index_vs[comps[*V].qtd_vs];
				comps[v2].index_vs[comps[v2].qtd_vs++] = aux;
				visitados[aux] = v2;
				comps[*V].index_vs[comps[*V].qtd_vs] = -1;
			}
		}
	}

}


void add_edges_of_label(int* visitados, Labels label, int *unique_edges, Problem* P, int(*m_adj)[P->V]){
	if(P->L == 0){
		P->edges_per_label = (Labels*) malloc((P->L+1) * sizeof(Labels));
	}else{
		P->edges_per_label = (Labels*) realloc(P->edges_per_label, (P->L+1) * sizeof(Labels));
	}
	P->edges_per_label[P->L].qtd_edges = 0;
	P->edges_per_label[P->L].es = NULL;
	P->edges_per_label[P->L].qtd_edge_per_k = malloc((P->k + 1)* sizeof(int));
	P->edges_per_label[P->L].my_value = label.my_value;
	int v1, v2, aux;
	for(int i = 0; i < label.qtd_edges; i++){
		v1 = label.es[i].v1;
		v2 = label.es[i].v2;
		if(visitados[v1] != visitados[v2]){
			v1 = visitados[v1];
			v2 = visitados[v2];
			if(m_adj[v1][v2] == -1){
				*unique_edges = *unique_edges + 1;
			}
			if(m_adj[v1][v2] != P->L){
				m_adj[v1][v2] = P->L;
				P->E++;
				add_edge(P, v1, v2, P->L);
			}
		}
	}
	for(int i = 0; i < P->k+1; i++)
		P->edges_per_label[P->L].qtd_edge_per_k[i] = P->edges_per_label[P->L].qtd_edges;
	P->L++;
}





void add_color_sub(VNS *vns, int r, int l_in_sub){
	vns->bool_labels[r] = vns->iteracao;
	vns->labels_sub[l_in_sub] = r; 
	// if(vns->ponderamento[r] > 1){
	// 	vns->ponderamento[r] -= 1;
	// 	vns->total_pond -= 1;
	// }
}

void create_vertices_and_matrix(int V, Problem *P, int(*m_adj)[P->V]){
	int i;

	P->edges_per_vertice = (Vertices*) malloc(V * sizeof(Vertices));
	for(i = 0; i < V; i++){
		P->edges_per_vertice[i].qtd_adj = 0;
		P->edges_per_vertice[i].es = NULL;
		for(int j = 0; j < V; j++){
			m_adj[i][j] = -1;
		}
	}
}

Problem 
gera_sub_h(VNS* vns, Problem G){	
	int visitados[G.V];
	COMPS *comps = init_comps(G.V, visitados);
	


	for(int v = 0; v < G.V; v++){
		visitados[v] = v;
	}

	int i, r, j, aux, V = G.V;


	int retake_total = vns->total_pond;
	for(i = 0; i < vns->qtd_fixed; i++){
		r = (rand() % retake_total) + 1;
		aux = 0;
		j = 0;
		while(aux < r){
			if(vns->bool_solution[j] != vns->iteracao){
				aux += vns->ponderamento[vns->solution_now[j]];
			}
			j++;
		}
		r = j-1;
		retake_total -= vns->ponderamento[vns->solution_now[r]];
		vns->bool_solution[r] = vns->iteracao; 
		r = vns->solution_now[r];
		vns->labels_fixed[i] = r;
		vns->bool_labels[r] = vns->iteracao;
		if(vns->ponderamento[r] > 1){
			vns->ponderamento[r] -= 1;
			vns->total_pond -= 1;
		}
		fix_labels(visitados, comps, G.edges_per_label[r], &V);
	}

	finish_comps(comps, G.V);



	Problem sub = start_problem();
	sub.V = V;
	sub.k = G.k - vns->qtd_fixed;
	int m_adj[V][V];
	create_vertices_and_matrix(V, &sub, m_adj);

	int unique_edges = 0;

	while(sub.L < floor((G.k - vns->qtd_fixed)*FRAC_SUB_FICA)){
		r = rand() % G.k;
		while(vns->bool_solution[r] == vns->iteracao){
			r = rand() % G.k;
		}
		vns->bool_solution[r] = vns->iteracao; 
		r = vns->solution_now[r];
		add_color_sub(vns, r, sub.L);
		add_edges_of_label(visitados, G.edges_per_label[r], &unique_edges, &sub, m_adj);
	}

	for(i = 0; i < G.k; i++){
		vns->bool_labels[vns->solution_now[i]] = vns->iteracao;
	}

	aux = 0;
	int ant_bug = 0;

	while(aux + G.k< FRAC_CORES_SUB_SAI *G.L){
		ant_bug = 0;

		while(vns->bool_labels[r] == vns->iteracao){
			r = rand() % G.L;
			ant_bug ++;
			if(ant_bug > G.L) break;
		}
		if(ant_bug > G.L) break;
		vns->bool_labels[r] = vns->iteracao;
		aux++;
	}

	for(i = 0; i < G.L; i++){
		if(vns->bool_labels[i] < vns->iteracao){
			add_color_sub(vns, i, sub.L);
			add_edges_of_label(visitados, G.edges_per_label[i], &unique_edges, &sub, m_adj);
		}
	}


	return(sub);
	
}





void count_arestas(int* visitados, Labels label, int V, int (*matrix)[V], int k_atual, int flag){
	int v1, v2, aux = 0;
	for(int i = 0; i < label.qtd_edges; i++){
		v1 = label.es[i].v1;
		v2 = label.es[i].v2;
		if(visitados[v1] != visitados[v2]){
			v1 = visitados[v1];
			v2 = visitados[v2];

			if(matrix[v1][v2] != flag){
				matrix[v1][v2] = flag;
				matrix[v2][v1] = flag;
				aux++; 
			}
		}
	}

	label.qtd_edge_per_k[k_atual] = aux;

}





void fix_labels2(int* visitados, COMPS* comps, Labels label, int *V, int(*edges)[2], int *prox_edge, int*edges_added, int k_atual){
	int v1, v2, aux;
	for(int i = 0; i < label.qtd_edges; i++){
		v1 = label.es[i].v1;
		v2 = label.es[i].v2;
		if(visitados[v1] != visitados[v2]){
			if(comps[visitados[v1]].qtd_vs > comps[visitados[v2]].qtd_vs){
				v2 = visitados[v2];
				v1 = visitados[v1];
			}else{
				aux = visitados[v2];
				v2 = visitados[v1];
				v1 = aux;
			}
			edges[*prox_edge][0] = v1;
			edges[*prox_edge][1] = v2;
			*prox_edge+=1;
			edges_added[k_atual+1]++;
			aux = comps[v2].qtd_vs;
			while(aux --){
				visitados[comps[v2].index_vs[aux]] = v1;
				comps[v1].index_vs[comps[v1].qtd_vs++] = comps[v2].index_vs[aux];
			}
			*V -= 1;
		}
	}

}


void defix_labels2(int k_atual, COMPS* comps, int* visitados, int(*edges)[2], int *last_edge, int* edges_added){
	int v1, v2, aux;
	while(edges_added[k_atual]){
		*last_edge-=1;
		v1 = edges[*last_edge][0];
		v2 = edges[*last_edge][1];
		for(aux = 0; aux < comps[v2].qtd_vs; aux++){
			comps[v1].qtd_vs--;
			visitados[comps[v1].index_vs[comps[v1].qtd_vs]] = v2;
			comps[v1].index_vs[comps[v1].qtd_vs] = -1;
		}
		edges[*last_edge][0] = -1;
		edges[*last_edge][1] = -1;
		edges_added[k_atual]--;
	}

}


int solve(Problem G, int(*edges)[2], int (*matrix)[G.V], int edge_atual, int *flag, COMPS* comps, int* visitados, int* best_value, int value_atual, int k_atual, int*solution_atual, int*best_solution, int*edges_added, int * solution_value){
	int v_fix, aux, i, j, l;
	

	for(l = solution_atual[k_atual] + 1; l < G.L; l++){

		v_fix = value_atual;
		fix_labels2(visitados, comps, G.edges_per_label[l], &v_fix, edges, &edge_atual, edges_added, k_atual);


		solution_atual[k_atual + 1] = l;
		solution_value[k_atual + 1] = G.edges_per_label[l].my_value;


		if(v_fix < *best_value){
			for(i = 1; i < k_atual + 2; i++){
				best_solution[i - 1] = solution_value[i];
				//printf("%d ", best_solution[i-1]);
			}
			*best_value = v_fix;
		}




		if(k_atual < G.k-1){
			for(i = l+1; i < G.L; i++){
				count_arestas(visitados, G.edges_per_label[i], G.V, matrix, k_atual + 1, *flag);
				*flag += 1;
			}


			mergeSort(G.edges_per_label, l + 1, G.L - 1, k_atual + 1);

		
			aux = 0;
			for(i = l + 1; i < l + 1 + G.k - k_atual && i < G.L; i++){
				aux += G.edges_per_label[i].qtd_edge_per_k[k_atual+1];
			}

			if(v_fix - aux < *best_value){
				solve(G, edges, matrix, edge_atual, flag, comps, visitados, best_value, v_fix, k_atual + 1, solution_atual, best_solution, edges_added, solution_value);
			}

		}

		defix_labels2(k_atual + 1, comps, visitados, edges, &edge_atual, edges_added);
	

		
	}


}

int generate_model_and_solve(Problem G, int *optimal){
	int best_value = G.V;
	int l, i, j;
	int solution_atual[G.k+1];
	int solution_value[G.k+1];
	int edges_added[G.k+1];
	int edges[G.V][2];
	int edge_atual = 0;
	int matrix[G.V][G.V];


	for(i = 0; i < G.V; i++){
		for(j = 0; j < G.V; j++){
			matrix[i][j] = -1;
		}
		edges[i][0] = -1;
		edges[i][1] = -1;
	}
	for(j = 0; j < G.k +1; j++){
		edges_added[j] = 0;
		solution_atual[j] = j;
		solution_value[j] = G.edges_per_label[j].my_value;
	}
	solution_atual[0] = -1;

	int flag = 0;
	int visitados[G.V];
	COMPS *comps = init_comps(G.V, visitados);
	mergeSort(G.edges_per_label, 0, G.L - 1, 1);
	solve(G, edges,  matrix, edge_atual, &flag, comps, visitados, &best_value, G.V, 0,  solution_atual, optimal, edges_added, solution_value);
	finish_comps(comps, G.V);
	return best_value;
}



void from_best_to_vns(VNS *vns, int k){
	vns->total_pond = 0;
	for(int i = 0; i < k; i++){
		vns->total_pond += vns->ponderamento[vns->best_know[i]];
		vns->solution_now[i] = vns->best_know[i];
	}
}

void 
from_solver_to_vns(VNS *vns, int *solution, int flag){
	int aux = vns->qtd_fixed;
	vns->total_pond = 0;
	for(int i = 0; i < vns->qtd_fixed; i++){
		vns->solution_now[i] = vns->labels_fixed[i];
		vns->total_pond +=  vns->ponderamento[vns->labels_fixed[i]];
		if(flag)
			vns->best_know[i] = vns->labels_fixed[i];
	}
	for(int i = 0; i < vns->qtd_labels_sub; i++){
		vns->solution_now[aux] = solution[i];
		vns->total_pond +=  vns->ponderamento[solution[i]];
		if(flag)
			vns->best_know[aux] = solution[i];
		aux++;

	}
}




int main(int argc, char **argv){

	int result = 0;

	FILE *saida;
	
	
	VNS vns;
	int error = 0;

	int it = 0;

	Problem G = start_problem();
	
	
	Problem sub = start_problem();

	error = read_file(&G, argv[1]);
	if(error) exit(1);

	FRAC_CORES_SUB_SAI = 1. - ((float) G.L / G.V);

	

	if(0.2 > FRAC_CORES_SUB_SAI){
		FRAC_CORES_SUB_SAI = 0.2;
	}


	
	
	
	srand(atoi(argv[2]));
	int mult = atoi(argv[3]);
	MAX_IN_SOLUTION = atoi(argv[4]); //1017
	FRAC_FIX = atof(argv[5]); 
	FRAC_SUB_FICA = atof(argv[6]);
	MAX_POND = atoi(argv[7]);
	int flaaag = atoi(argv[8]);
	
	MAX_ITERATION = G.k * mult;

	int sol_to_avaliation[G.L];
	vns = start_alg(G.V, G.L, G.k);

	int solution_solver[G.L];
	for(int i = 0; i < G.L; i++){
		solution_solver[i] = 0;
	}
	result = vns.best;
	if(vns.best > 1 && (time(NULL) - vns.time) < 300){
		while(it < MAX_ITERATION && (time(NULL) - vns.time) < 300){
			vns.iteracao++;

			sub = gera_sub_h(&vns,  G);

			//printa_infos(vns, sub.V, sub.E, sub.L, sub.k);
			result = generate_model_and_solve(sub, solution_solver);

			if(result <= vns.solution_value && result > 0.1){
				vns.solution_value = result;
				if(vns.solution_value < vns.best){	
					from_solver_to_vns(&vns, solution_solver, 1);
					vns.best = vns.solution_value;

					if(vns.best == 1){
						break;
					}
					//vns.max_in_solution = MAX_IN_SOLUTION;

				}else{
					from_solver_to_vns(&vns, solution_solver, 0);	
				}
				vns.it_in_solution = 0;
				
				
			}else{
				vns.it_in_solution ++;
				if(vns.it_in_solution % vns.max_in_solution == 0){
					//vns.max_in_solution += MAX_IN_SOLUTION;
					vns.solution_value = result;
					from_solver_to_vns(&vns, solution_solver, 0);
				}
				// else if(((vns.it_in_solution + 1) % (MAX_IN_SOLUTION * 5) == 0)){
				// 	from_best_to_vns(&vns, G.k);
				// 	vns.it_in_solution = 0;
				// }
			}
			if(it % 50 == 0){
				printf("\n%d %d", vns.solution_value, vns.it_in_solution);
				printf("\n%s\nbest: %d\nit: %d \ntempo: %ld \n\n", argv[1], vns.best, vns.iteracao, time(NULL) - vns.time);
			}

			it+=1;		
			finish_problem(&sub);
		}

	

		
		
	}

	for(int i = 0; i < G.L; i++)
		sol_to_avaliation[i] = 0;
	for(int i = 0; i < G.k; i ++)
		sol_to_avaliation[vns.best_know[i]] = 1;
	vns.best = avalia_solution(sol_to_avaliation, G);
	// if(aux != vns.best ){
	// 	printf("\n%d %d\n", aux, vns.best);

	// 	finish_problem(&G); 
	// 	finish_vns(&vns);
	// 	exit(1);
	// }


	if(flaaag){
		char file_saida[50];
		snprintf(file_saida, sizeof(file_saida), "saidas/saida_%d_%d_%0.2f_%0.2f_%d", mult, MAX_IN_SOLUTION, FRAC_FIX, FRAC_SUB_FICA, MAX_POND);
		saida = fopen(file_saida, "a");	
		fprintf(saida, "%s:%ld:%d:1\n", argv[1], time(NULL) - vns.time, vns.best);	
		fclose(saida);
	}
	printf("\nBest %d", vns.best);
	finish_problem(&G); 
	finish_vns(&vns);

	exit(0);
}
