#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "gurobi_c.h"
#include <string.h>
#include <math.h>

// export GUROBI_HOME="gurobi912/linux64"
// export PATH="${PATH}:${GUROBI_HOME}/bin"
// export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"

// gcc -o p v_restriction_extra.c -I gurobi912/linux64/include/ -L gurobi912/linux64/lib/ -lgurobi91 -lm


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
	int qtd_edges;
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
	Labels *edges_per_label;
	Vertices *edges_per_vertice;
} Problem;

typedef struct COMPS{
	int qtd_vs;
	int *index_vs;
}COMPS;




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
	int it_find_best;
	GRBenv *env;
	double best;
	double solution_value;
	long time;
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

	G->edges_per_label = (Labels*) malloc(G->L*sizeof(Labels));
	if(G->edges_per_label == NULL){
		printf("Out of memory");
		exit(1);
	}

	for (int i = 0; i < G->L; i++){
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
	GRBfreeenv(vns->env);
}

void 
finish_problem(Problem *G){
	if(G->edges_per_label != NULL){
		for (int l = 0; l < G->L; l ++){
			if(G->edges_per_label[l].es != NULL)
				free(G->edges_per_label[l].es);
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

float 
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
	fprintf(saida, "\n%.0lf %.0lf %ld %d %d %d", vns.best, vns.solution_value, vns.time, vns.iteracao, vns.it_in_solution, vns.total_pond);
	fclose(saida);
}


void 
errorr(GRBenv *env){
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(1);
}








COMPS *init_comps(int V){
	COMPS* comps = (COMPS*) malloc(V * sizeof(COMPS));;
	for(int v = 0; v < V; v++){
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







void 
generate_first_solution(VNS *vns, Problem G){
	int min = G.V, id_min = -1, i = 0, j =0;
	
	int solution[G.L];

	for(i = 0; i < G.L; i++){
		solution[i] = 0;
	}

	vns->best = G.V;
	double result;

	for(i = 0; i < G.k; i++ ){
		result = vns->best;
		for(j = 0; j < G.L; j++){
			if(solution[j] == 0){
				solution[j] = 1;
				result = avalia_solution(solution, G);
				if(min > result){
					id_min = j;
					min = result;
				}
				solution[j] = 0;
			}
		}
		solution[id_min] = 1;
		vns->bool_labels[id_min] = vns->iteracao;
		vns->solution_now[i] = id_min;
		vns->bool_solution[i] = -1;
		vns->best_know[i] = id_min;
		vns->best = min;
	}

}


void printa_comps(COMPS *comps, int V){
	for(int i = 0; i < V; i++){
		printf("\n %d: ", comps[i].qtd_vs);
		for(int j = 0; j < V; j++){
			printf("%d ", comps[i].index_vs[j]);
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


Problem 
gera_sub_h(VNS* vns, Problem G){	
	COMPS *comps = init_comps(G.V);
	
	int visitados[G.V];


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



	Problem sub;
	sub.V = V;
	sub.E = 0;
	sub.L = 0;
	sub.k = G.k - vns->qtd_fixed;



	int m_adj[V][V];
	sub.edges_per_vertice = (Vertices*) malloc(V * sizeof(Vertices));
	for(i = 0; i < V; i++){
		sub.edges_per_vertice[i].qtd_adj = 0;
		sub.edges_per_vertice[i].es = NULL;
		for(int j = 0; j < V; j++){
			m_adj[i][j] = -1;
		}
	}

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
	while(aux + G.k < FRAC_CORES_SUB_SAI *G.L){
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

	vns->qtd_labels_sub = sub.L;
	aux = sub.L;
	sub.L += sub.V - 1;
	sub.edges_per_label = (Labels*) realloc(sub.edges_per_label, sub.L * sizeof(Labels));
	for(int v = 0; v < sub.V - 1; v++)
	{		
		sub.edges_per_label[aux + v].qtd_edges = 0;
		sub.edges_per_label[aux + v].es = NULL;
		add_edge(&sub, sub.V - 1, v, aux + v);
		sub.E += 1;
	}
	return(sub);
	
}



int 
DFS_gera_sub(int v, int* visitados,  double* sorted_labels, Problem *sub, Problem G, int(*bool_adjs_labels)[G.L]){
	visitados[v] = sub->V - 1;
	int v1 = sub->V-1;
	int v2, l;
	Edges e;

	for(int i = 0; i < G.edges_per_vertice[v].qtd_adj; i++){
		e = G.edges_per_vertice[v].es[i];
		if(e.v1 == v){
			v2 = e.v2;
			l = e.l;
		}else{
			v2 = e.v1;
			l = e.l;
		}
		if(sorted_labels[l] == -1 && visitados[v2] == -1){
			DFS_gera_sub(v2, visitados, sorted_labels, sub, G, bool_adjs_labels);
		}else if(visitados[v2] != v1 && visitados[v2] != -1){
			v2 = visitados[v2];
			l = sorted_labels[l];
			if(bool_adjs_labels[v2][l] != v1){
				bool_adjs_labels[v2][l] = v1;
				add_edge(sub, v1, v2, l);
				sub->E++;
			}
		}
	}
}

Problem 
gera_sub_s(double* solution, Problem G, double teste, int*ordem){
	int aux = 0;
	for(int l = 0; l < G.L; l++){
		if(solution[l] > teste)
			solution[l] = -1;
		else{
			ordem[aux] = l;
			solution[l] = aux++;
		}
	}


	Problem sub;
	sub.V = 0;
	sub.E = 0;
	sub.L = aux;
	sub.k = 0;

	sub.edges_per_label = (Labels*) malloc(sub.L* sizeof(Labels));
	sub.edges_per_vertice = NULL;
	for(int l = 0; l < sub.L; l++){
		sub.edges_per_label[l].qtd_edges = 0;
		sub.edges_per_label[l].es = NULL;
	}
	

	int visitados[G.V];
	int bool_adjs_labels[G.V][sub.L];
	for(int v = 0; v < G.V; v++){
		visitados[v] = -1;
		for(int l = 0; l < sub.L; l++)
			bool_adjs_labels[v][l] = -1;
	}

	
	for (int v = 0; v < G.V; v++){
		if(visitados[v] == -1){
			sub.V += 1;
			if(sub.V == 1)
				sub.edges_per_vertice = (Vertices*) malloc(sub.V * sizeof(Vertices));
			else
				sub.edges_per_vertice = (Vertices*) realloc(sub.edges_per_vertice, sub.V * sizeof(Vertices));
			sub.edges_per_vertice[sub.V-1].qtd_adj = 0;
			sub.edges_per_vertice[sub.V-1].es = NULL;
			DFS_gera_sub(v, visitados, solution, &sub, G, bool_adjs_labels);
		}
	}
	
	return(sub);
}


Problem 
gera_sub_f(double* bool_labels, Problem G, double teste){
	COMPS *comps = init_comps(G.V);
	
	int visitados[G.V];



	for(int v = 0; v < G.V; v++){
		visitados[v] = v;
	}

	int aux = 0, V = G.V;

	for(int i = 0; i < G.L; i++){
		if(bool_labels[i] > teste){
			fix_labels(visitados, comps, G.edges_per_label[i], &V);
		}
	}

	finish_comps(comps, V);

	Problem sub;
	sub.V = V;
	sub.E = 0;
	sub.L = 0;
	sub.k = 0;
	sub.edges_per_label = NULL;


	int m_adj[V][V];
	sub.edges_per_vertice = (Vertices*) malloc(V * sizeof(Vertices));
	for(int i = 0; i < V; i++){
		sub.edges_per_vertice[i].qtd_adj = 0;
		sub.edges_per_vertice[i].es = NULL;
		for(int j = 0; j < V; j++){
			m_adj[i][j] = -1;
		}
	}

	int unique_edges = 0;
	for(int i = 0; i < G.L; i++){
		if(round(bool_labels[i]) < teste){
			bool_labels[sub.L] = i;
			add_edges_of_label(visitados, G.edges_per_label[i], &unique_edges, &sub, m_adj);
		}
	}



	return(sub);
}




int gera_cortes(Problem sub,  void* cbdata, int *ordem, int type, int L){
	int len_restrict = 0;
	int *ind = NULL;
	double *val = NULL;
	int l = 0;
	int error = 0;
	int present_in_cort[L];
	for(int i = 0; i < L; i++){
		present_in_cort[i] = -1;
	}

	int j;


	if(sub.V > 1){
		for(int v = 0; v < sub.V; v++){
			len_restrict = 0;
			for( j = 0; j < sub.edges_per_vertice[v].qtd_adj; j++){
				l = sub.edges_per_vertice[v].es[j].l;
				l = ordem[l];
				if(present_in_cort[l] != v){
					present_in_cort[l] = v;
					len_restrict++;
					if(len_restrict == 1){
						ind = (int*) malloc(len_restrict*sizeof(int));
						val = (double*) malloc(len_restrict*sizeof(double));
					}else{
						ind = (int*) realloc(ind, len_restrict*sizeof(int));
						val = (double*) realloc(val, len_restrict*sizeof(double));
					}
					ind[len_restrict-1] = l;
					val[len_restrict-1] = 1;
				}
			}
			if(len_restrict){
				if(type)
					error = GRBcblazy(cbdata, len_restrict, ind, val, GRB_GREATER_EQUAL, 1.);
				else
					error = GRBcbcut(cbdata, len_restrict, ind, val, GRB_GREATER_EQUAL, 1.);
				free(ind);
				free(val);
			}
		}
	}
	return(error);
}



int __stdcall 
calback_func(GRBmodel *model, void *cbdata, int where, void *infos){
	Problem *master_h = (Problem*) infos;

	Problem sub = start_problem();
	double sol[master_h->L];
	int ordem[master_h->L];
	for(int i = 0; i < master_h->L; i++){
		sol[i] = 0;
		ordem[i] = -1;
	}

	int error = 0;

	if (where == GRB_CB_MIPSOL) {
	    GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol);
		sub = gera_sub_s(sol, *master_h, 0.1, ordem);
		error = gera_cortes(sub, cbdata, ordem, 1, master_h->L);
		finish_problem(&sub);
	}else if (where == GRB_CB_MIPNODE) {
		GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol);
		sub = gera_sub_s(sol, *master_h, 0.0, ordem);
		error = gera_cortes(sub, cbdata, ordem, 0, master_h->L);
		finish_problem(&sub);
	}
	return(error);
}



double generate_model_and_solve(int qtd_labels_sub, Problem sub_h, double*solution, GRBenv *env, int best){
	char name[100];


	int error;
	GRBmodel *model = NULL;
	int ind[qtd_labels_sub];
	double val[qtd_labels_sub];
	double val2[sub_h.L];
	int ind2[sub_h.L];
	double val_extra[sub_h.L - qtd_labels_sub];
	int ind_extra[sub_h.L - qtd_labels_sub];
	int solcount;
	int r;
  	double result = sub_h.V;
	

	
	error = GRBnewmodel(env, &model, "klsf", 0, NULL, NULL, NULL, NULL, NULL);
	if (error) errorr(env);
	

	for(int l = 0; l < sub_h.L; l++){
		if(l < qtd_labels_sub){
			solution[l] = 0;
			ind[l] = l;
			val[l] = 1.0;
			sprintf(name, "x_%d", l);
			error = GRBaddvar(model, 0, NULL, NULL, 0.0 , 0.0, 1.0, GRB_BINARY, name);
		}else{
			sprintf(name, "xx_%d", l);	
			error = GRBaddvar(model, 0, NULL, NULL, 1.0 , 0.0, 1.0, GRB_BINARY, name);
			if(best){
				val_extra[l - qtd_labels_sub] = 1;
				ind_extra[l - qtd_labels_sub] = l;
			}
		}
		val2[l] = sub_h.edges_per_label[l].qtd_edges;
		ind2[l] = l;
		if (error) errorr(env);
	}

	sprintf(name, "extra");
	if(best){
		error = GRBaddconstr(model, sub_h.L - qtd_labels_sub, ind_extra, val_extra, GRB_LESS_EQUAL, best, name);
		if (error) errorr(env);
	}
	sprintf(name, "cores");
	error = GRBaddconstr(model, qtd_labels_sub, ind, val, GRB_LESS_EQUAL, sub_h.k, name);
	if (error) errorr(env);
  	sprintf(name, "v-1");
	error = GRBaddconstr(model, sub_h.L, ind2, val2, GRB_GREATER_EQUAL, sub_h.V-1, name);
	if (error) errorr(env);
  	
  	//error = GRBwrite(model, "model.rlp.gz");
  	//if(error) errorr(env);	
	error = GRBsetcallbackfunc(model, calback_func, (void *) &sub_h);
	if (error) errorr(env);
	error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LAZYCONSTRAINTS, 1);
	if (error) errorr(env);
	error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
	if (error) errorr(env);
	error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", 3.0);
	if (error) errorr(env);
	error = GRBsetintparam(GRBgetenv(model), "Threads", 1);
	if (error) errorr(env);
	long int t_solver = time(NULL);
	error = GRBoptimize(model);
//	if(time(NULL) - t_solver > 2){
//		FRAC_SUB_FICA += 0.05;
//		MAX_IN_SOLUTION += 5;
//	}
	if (error) errorr(env);
	error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount);
  	if (error) errorr(env);
 	
 	if(solcount){
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &result);
  		if (error) errorr(env);
  		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, qtd_labels_sub, solution);
  		if (error) errorr(env);
  	}else{
  		result = sub_h.V;
  		for(int s = 0; s < sub_h.k; s++){
  			r = rand()%qtd_labels_sub;
  			while(solution[r] >= 1){
  				r = rand()%qtd_labels_sub;	
  			}
  			solution[r] = 1;
  		}

  	}



	GRBfreemodel(model);

	return(result);
}

void from_best_to_vns(VNS *vns, int k){
	vns->total_pond = 0;
	for(int i = 0; i < k; i++){
		vns->total_pond += vns->ponderamento[vns->best_know[i]];
		vns->solution_now[i] = vns->best_know[i];
	}
}


void 
from_solver_to_vns(VNS *vns, double *solution, int flag){
	int aux = vns->qtd_fixed;
	vns->total_pond = 0;
	for(int i = 0; i < vns->qtd_fixed; i++){
		vns->solution_now[i] = vns->labels_fixed[i];
		vns->total_pond +=  vns->ponderamento[vns->labels_fixed[i]];
		if(flag)
			vns->best_know[i] = vns->labels_fixed[i];
	}
	for(int i = 0; i < vns->qtd_labels_sub; i++){
		if(round(solution[i]) > 0.9){
			vns->solution_now[aux] = vns->labels_sub[i];
			vns->total_pond +=  vns->ponderamento[vns->labels_sub[i]];
			if(flag)
				vns->best_know[aux] = vns->labels_sub[i];
			aux++;
		}

	}
}


VNS 
start_alg(Problem G){

	VNS vns;
	vns.it_find_best = 0;
	int error = GRBloadenv(&vns.env, "solver.log");
	if (error) errorr(vns.env);
	vns.max_in_solution = MAX_IN_SOLUTION;
	vns.iteracao = 0;
	vns.ponderamento = (int*) malloc(G.L * sizeof(int));
	vns.best_know = (int*) malloc(G.k*sizeof(int));
	vns.bool_labels = (int*) malloc(G.L*sizeof(int));
	vns.qtd_labels_sub = 0;
	vns.qtd_fixed = (int) floor(G.k * FRAC_FIX);
	if(G.k == vns.qtd_fixed){
		vns.qtd_fixed --;
	}
	vns.labels_fixed = (int*) malloc(vns.qtd_fixed * sizeof(int));
	vns.labels_sub = (int*) malloc(G.L * sizeof(int));
	vns.solution_now = (int*) malloc(G.k*sizeof(int));
	vns.bool_solution = (int*) malloc(G.k*sizeof(int));
	vns.best = G.V;
	vns.solution_value = G.V;
	vns.time = time(NULL);
	vns.it_in_solution = 0;


	vns.total_pond = 0;
	for(int i = 0; i < G.L ; i++){
		vns.labels_sub[i] = -1;
	}

	for(int i = 0; i < vns.qtd_fixed; i++){
		vns.labels_fixed[i] = 0;
	}

	for(int i = 0; i < G.L; i++){
		vns.ponderamento[i] = MAX_POND;
		vns.bool_labels[i] = -1;
	}
	vns.total_pond = G.k * MAX_POND;

	int cont_l = 0, r;

	generate_first_solution(&vns, G);
	
	return(vns);
}




int main(int argc, char **argv){

	double result = 0;

	FILE *saida;
	
	
	VNS vns;
	int error = 0;

	int it = 0;

	Problem G = start_problem();
	
	
	Problem sub = start_problem();

	error = read_file(&G, argv[1]);
	if(error) exit(1);

	FRAC_CORES_SUB_SAI = 1. - ((double) G.L / G.V);

	

	if(0.2 > FRAC_CORES_SUB_SAI){
		FRAC_CORES_SUB_SAI = 0.2;
	}


	
	
	
	srand(atoi(argv[2]));
	int mult = atoi(argv[3]);
	float mult2 = atof(argv[4]); //1017
	FRAC_FIX = 0.8;
	MAX_POND = atoi(argv[5]);
	int flaaag = atoi(argv[6]);
	MAX_IN_SOLUTION = floor(G.k * mult2);
	MAX_ITERATION = G.k * mult;

	int sol_to_avaliation[G.L];
	vns = start_alg(G);

	double solution_solver[G.L];
	for(int i = 0; i < G.L; i++){
		solution_solver[i] = 0;
	}
	result = vns.best;
	if(vns.best > 1 && (time(NULL) - vns.time) < 300){
		while(it < MAX_ITERATION && (time(NULL) - vns.time) < 300){
			vns.iteracao++;
			
			if(vns.it_in_solution % vns.max_in_solution == 0)
				FRAC_SUB_FICA = 0.1;
			else
				FRAC_SUB_FICA = 0.8;

			
			sub = gera_sub_h(&vns,  G);

			//printa_infos(vns, sub.V, sub.E, sub.L, sub.k);
			//printf("\n");



			if(vns.it_in_solution % vns.max_in_solution == 0)
				result = generate_model_and_solve(vns.qtd_labels_sub, sub, solution_solver, vns.env, 0);
			else
				result = generate_model_and_solve(vns.qtd_labels_sub, sub, solution_solver, vns.env, vns.solution_value);

			result = round(result) + 1; 

			if(result <= vns.solution_value && result > 0.1){
				vns.solution_value = round(result);
				if(vns.solution_value < vns.best){	
					from_solver_to_vns(&vns, solution_solver, 1);
					vns.best = vns.solution_value;
					vns.it_find_best = it;
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
				printf("\n%.0lf %d", vns.solution_value, vns.it_in_solution);
				printf("\n%s\nbest: %.0f\nit: %d \ntempo: %ld \n\n", argv[1], vns.best, vns.iteracao, time(NULL) - vns.time);
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
	// 	printf("\n%d %lf\n", aux, vns.best);

	// 	finish_problem(&G); 
	// 	finish_vns(&vns);
	// 	exit(1);
	// }


	if(flaaag){
		char file_saida[50];
		snprintf(file_saida, sizeof(file_saida), "saidas/s_%d_%.2lf_%.2f_%d", mult, mult2,  FRAC_FIX, MAX_POND);
		saida = fopen(file_saida, "a");	
		fprintf(saida, "%s:%ld:%.0lf:%d\n", argv[1], time(NULL) - vns.time, vns.best, vns.it_find_best);	
		fclose(saida);
	}
	printf("\nBest %f", vns.best);
	finish_problem(&G); 
	finish_vns(&vns);

	exit(0);
}



