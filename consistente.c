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

//gcc -o p consistente.c -I gurobi912/linux64/include/ -L gurobi912/linux64/lib/ -lgurobi91 -lm


#define MAX_ITERATION 50
#define TOTAL_SUB 20.0
#define FRAC_FIX 0.55

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


typedef struct VNS{
	int *best_know;
	int *solution_now;
	int *ponderamento;
	int *list_labels_sub;
	int *mapped_label;
	int iteracao;
	int qtd_labels_sub;
	int qtd_fixed;
	int total_pond;
}VNS;


void 
input(){
	int x;
	printf("\n");
	scanf("%d", &x);
}

void 
printa_infos(VNS vns, int V, int E, int L, int k){
	//printf("\nlist_labels_sub: ");
	//for(int i = 0; i < vns.qtd_labels_sub; i++)
	//	printf("%d ", vns.list_labels_sub[i]);
	// printf("\nsolution_now: ");
	// for(int i = 0; i < k; i++)
	// 	printf("%d ", vns.solution_now[i]);
	// printf("\nmapped_label: ");
	// for(int i = 0; i < L; i++)
	// 	printf("%d ", vns.mapped_label[i]);
	// printf("\nbest_know: ");
	// for(int i = 0; i < k; i++)
	// 	printf("%d ", vns.best_know[i]);
	printf("\nqtd_labels_sub: %d", vns.qtd_labels_sub);

}


void 
printa_problem(Problem P){
	printf("\nV: %d, E: %d L: %d, k: %d", P.V, P.E, P.L, P.k);
	// for (int l = 0; l < P.L; l++){
	// 	printf("\n%d: ", l);
	// 	Edges e;
	// 	for(int i = 0; i < P.edges_per_label[l].qtd_edges; i++){
	// 		e = P.edges_per_label[l].es[i];
	// 		printf("(%d, %d)", e.v1, e.v2);
	// 	}
	// }
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

	printf("%d %d %d %d \n", G->V, G->E, G->L, G->k);

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
	if(vns->list_labels_sub != NULL)
		free(vns->list_labels_sub);
	if(vns->solution_now != NULL)
		free(vns->solution_now);
	if(vns->mapped_label != NULL)
		free(vns->mapped_label);
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

int 
avalia_solution(int *solution, Problem G){
	int componente = 0.;

	int visitados[G.V];

	for(int i = 0; i < G.V; i++)
		visitados[i] = -1;

	for(int v = 0; v < G.V; v++){
		if(visitados[v] == -1){
			componente =  componente + 1;
			DFS_avalia(v, componente, solution, visitados, G);
		}
	}

	return(componente);
}

VNS 
start_alg(int V, int E, int L, int k){
	VNS vns;
	vns.iteracao = 0;
	vns.ponderamento = (int*) malloc(L * sizeof(int));
	vns.best_know = (int*) malloc(k*sizeof(int));
	vns.mapped_label = (int*) malloc(L*sizeof(int));
	vns.qtd_labels_sub = (int) ceil(k * TOTAL_SUB);
	vns.qtd_fixed = (int) ceil(k * FRAC_FIX);
	vns.list_labels_sub = (int*) malloc(vns.qtd_labels_sub * sizeof(int));
	vns.solution_now = (int*) malloc(k*sizeof(int));

	for(int i = 0; i < k; i++){
		vns.solution_now[i] = -1;
		vns.best_know[i] = -1;
	}
	vns.total_pond = 0;
	for(int i = 0; i < vns.qtd_labels_sub; i++){
		vns.list_labels_sub[i] = -1;
	}
	for(int i = 0; i < L; i++){
		vns.ponderamento[i] = 99;
		vns.total_pond += 99;
		vns.mapped_label[i] = -1;
	}
	int cont_l = 0, r;
	while(cont_l < k){
		r = rand()%L;
		while(vns.mapped_label[r] != -1){
			r = rand()%L;
		}
		vns.solution_now[cont_l] = r;
		vns.best_know[cont_l] = r;
		cont_l += 1;
	}
	return(vns);
}






void 
pertubation_and_fix_values(VNS *vns, int k, int L){
	int r;
	int aux, j;
	vns->iteracao+=1;

	for(int i = 0; i < vns->qtd_labels_sub; i++){
		if (i < vns->qtd_fixed){
			r = rand() % k;
			while(vns->solution_now[r] == -1)
				r = rand() % k;
			aux = r;
			r = vns->solution_now[r];
			vns->solution_now[aux] = -1;
		}else{
			while(vns->mapped_label[r] >= 0){
//				r = rand() % L;
				r = (rand() % vns->total_pond) + 1;
				aux = 0;
				j = 0;
				while(aux < r){
					aux += vns->ponderamento[j];
					j+=1;					
				}
				r = j-1;
			}
		}
		
		vns->list_labels_sub[i] = r; 
		vns->mapped_label[r] = i; 
		if(vns->ponderamento[r] > 1){
			vns->ponderamento[r] -= 1;
			vns->total_pond -= 1;
		}
	}
}

void 
from_solver_to_vns(VNS *vns, double *solution){
	int aux = vns->qtd_fixed;
	for(int i = 0; i < vns->qtd_labels_sub; i++){
		if(i < vns->qtd_fixed){
			vns->solution_now[i] = vns->list_labels_sub[i];
			vns->best_know[i] = vns->list_labels_sub[i];
		}else{
			if(round(solution[i - vns->qtd_fixed]) > 0.9){
				vns->solution_now[aux] = vns->list_labels_sub[i];
				vns->best_know[aux] = vns->list_labels_sub[i];
				aux ++;
			}
		}
		vns->mapped_label[vns->list_labels_sub[i]] = -1;
		vns->list_labels_sub[i] = -1;
	}	
}


void 
random_solution_teste(VNS *vns, int k){
	int index, index_k = 0;
	int r;

	for(int i = 0; i < vns->qtd_labels_sub - k; i++){
		r = rand()%vns->qtd_labels_sub;
		r = vns->list_labels_sub[r];
		while(vns->mapped_label[r] >= 0){
			r = rand()%vns->qtd_labels_sub;
			r = vns->list_labels_sub[r];
		}
		vns->mapped_label[r] = vns->iteracao - 1;
	}
}


int 
DFS_gera_sub(int v, int* visitados,  int* mapped_label, int qtd_fixed, Problem *sub_h, Problem G, int(*bool_adjs_labels)[G.L], int*bool_adjs_edges,  int*qtd_unique_edges){
	visitados[v] = sub_h->V - 1;
	int v1 = sub_h->V-1;
	int v2, l;
	Edges e;

	for(int i = 0; i < G.edges_per_vertice[v].qtd_adj; i++){
		e = G.edges_per_vertice[v].es[i];
		if(mapped_label[e.l] != -1){
			if(e.v1 == v){
				v2 = e.v2;
				l = e.l;
			}else{
				v2 = e.v1;
				l = e.l;
			}
			if(mapped_label[l] < qtd_fixed && visitados[v2] == -1){
				DFS_gera_sub(v2, visitados, mapped_label, qtd_fixed, sub_h, G, bool_adjs_labels, bool_adjs_edges, qtd_unique_edges);
			}else if(visitados[v2] != v1 && visitados[v2] != -1){
				v2 = visitados[v2];
				l = mapped_label[l] - qtd_fixed;
				if(bool_adjs_labels[v2][l] != v1){
					if(bool_adjs_edges[v2] != v1){
						bool_adjs_edges[v2] = v1;
						qtd_unique_edges++;
					}
					bool_adjs_labels[v2][l] = v1;
					add_edge(sub_h, v1, v2, l);
					sub_h->E++;
				}
			}
		}
	}
}

Problem 
gera_sub_h(int* mapped_label, int qtd_labels_sub, int qtd_fixed, int flag_extend, Problem G){
	Problem sub_h;
	sub_h.V = 0;
	sub_h.E = 0;
	sub_h.L = qtd_labels_sub - qtd_fixed;
	sub_h.k = G.k - qtd_fixed;
	int qtd_unique_edges = 0;
	sub_h.edges_per_label = (Labels*) malloc(sub_h.L* sizeof(Labels));
	for(int l = 0; l < sub_h.L; l++){
		sub_h.edges_per_label[l].qtd_edges = 0;
		sub_h.edges_per_label[l].es = NULL;
	}
	sub_h.edges_per_vertice = NULL;

	int visitados[G.V];
	int bool_adjs_edges[G.V];
	int bool_adjs_labels[G.V][G.L];
	for(int v = 0; v < G.V; v++){
		visitados[v] = -1;
		bool_adjs_edges[v] = -1;
		for(int l = 0; l < sub_h.L; l++)
			bool_adjs_labels[v][l] = -1;

	}

	
	for (int v = 0; v < G.V; v++){
		if(visitados[v] == -1){
			sub_h.V += 1;
			if(sub_h.V == 1)
				sub_h.edges_per_vertice = (Vertices*) malloc(sub_h.V * sizeof(Vertices));
			else
				sub_h.edges_per_vertice = (Vertices*) realloc(sub_h.edges_per_vertice, sub_h.V * sizeof(Vertices));
			sub_h.edges_per_vertice[sub_h.V-1].qtd_adj = 0;
			sub_h.edges_per_vertice[sub_h.V-1].es = NULL;
			DFS_gera_sub(v, visitados, mapped_label, qtd_fixed, &sub_h, G, bool_adjs_labels, bool_adjs_edges, &qtd_unique_edges);
		}
	}
	if(flag_extend){
		int aux = sub_h.L;
		sub_h.L =  aux + sub_h.V;
		sub_h.edges_per_label = (Labels*) realloc(sub_h.edges_per_label, sub_h.L * sizeof(Labels));
		sub_h.V += 1;
		sub_h.edges_per_vertice = (Vertices*) realloc(sub_h.edges_per_vertice, sub_h.V * sizeof(Vertices));
		sub_h.edges_per_vertice[sub_h.V-1].qtd_adj = 0;
		sub_h.edges_per_vertice[sub_h.V-1].es = NULL;
		for(int v = 0; v < sub_h.V - 1; v++){
			sub_h.edges_per_label[aux + v].qtd_edges = 0;
			sub_h.edges_per_label[aux + v].es = NULL;
			add_edge(&sub_h, sub_h.V - 1, v, aux + v);
			sub_h.E += 1;
		}
	}

	return(sub_h);
}


int 
remove_labels(int* first_solution, int len_first, Problem *master){
	int aux = -1;
    Problem sub = start_problem();
    Problem max_problem;
   	int next_sol[master->L];
    int too_remove = -1;

    while(master->V >= 2){
    	max_problem = start_problem();
    	for(int i = 0; i < master->L; i++){
			next_sol[i] = i;
		}	
		too_remove = -1;
		for(int i = 0; i < master->L; i++){
			next_sol[i] = 0;
			sub = gera_sub_h(next_sol, master->L, 1, 0, *master);
			if(sub.V >= max_problem.V){
				finish_problem(&max_problem);
				max_problem = sub;
				too_remove = i;
			}else{
				finish_problem(&sub);
			}
			next_sol[i] = i + 1;
		}
		aux = -1;
		if(too_remove >= 0){
			for(int j = 0; j < len_first; j ++){
				if (first_solution[j] >= 0){
					aux++;
					if(aux == too_remove){
						first_solution[j] = -2;
						aux = j;	
					}
				}
			}
		}
		finish_problem(master);
		*master = max_problem;	
	}
	finish_problem(master);

	if(aux >= 0)
		first_solution[aux] = 1; 
}



int __stdcall 
calback_func(GRBmodel *model, void *cbdata, int where, void *infos){
	Problem *master_h = (Problem*) infos;

	Problem sub = start_problem();
	double sol[master_h->L];
	for(int i = 0; i < master_h->L; i++)
		sol[i] = 0;
	int labels_sub[master_h->L];
	int bool_vector_l[master_h->L];
	int aux = 0,  len_restrict = 0;
	int ordem[master_h->L];
	int *ind = NULL;
	double *val = NULL;
	int error = 0;
	Edges e;
	int l = 0;
	if (where == GRB_CB_MIPSOL) {
	    GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol);
		for(int i = 0; i < master_h->L; i++){
			ordem[i] = -2;
			if(sol[i] > 0.9){
				labels_sub[i] = -2;
			}else{
				ordem[aux] = i;
				labels_sub[i] = aux++;
			}
			bool_vector_l[i] = -1;
		}
		sub = gera_sub_h(labels_sub, aux, 0, 0, *master_h);

		if(sub.V > 1){
			for(int v = 0; v < sub.V; v++){
				for(int j = 0; j < sub.edges_per_vertice[v].qtd_adj; j++){
					e = sub.edges_per_vertice[v].es[j];
					l = ordem[e.l];
					if(bool_vector_l[l] != v){
						bool_vector_l[l] = v;
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
					error = GRBcblazy(cbdata, len_restrict, ind, val, GRB_GREATER_EQUAL, 1.);
					free(ind);
					ind = NULL;
					free(val);
					val = NULL;
					len_restrict = 0;
				}
			}
		}
		finish_problem(&sub);
	}else if (where == GRB_CB_MIPNODE) {
		GRBcbget(cbdata, where, GRB_CB_MIPSOL_SOL, sol);
		for(int i = 0; i < master_h->L; i++){
			ordem[i] = -2;
			if(sol[i] > 0){
				labels_sub[i] = -2;
			}else{
				ordem[aux] = i;
				labels_sub[i] = aux++;
			}
			bool_vector_l[i] = -1;
		}
		sub = gera_sub_h(labels_sub, aux, 0, 0, *master_h);
		if(sub.V > 1){
			for(int v = 0; v < sub.V; v++){
				len_restrict = 0;
				for(int j = 0; j < sub.edges_per_vertice[v].qtd_adj; j++){
					e = sub.edges_per_vertice[v].es[j];
					l = ordem[e.l];
					if(bool_vector_l[l] != v){
						bool_vector_l[l] = v;
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
					error = GRBcbcut(cbdata, len_restrict, ind, val, GRB_GREATER_EQUAL, 1.);
					free(ind);
					ind = NULL;
					free(val);
					val = NULL;
					len_restrict = 0;
				}
			}
		}
		finish_problem(&sub);
	}
	return(error);
}

void 
errorr(GRBenv *env){
	printf("ERROR: %s\n", GRBgeterrormsg(env));
	exit(1);
}

double* generate_model_and_solve(Problem sub_h, double *result, int V){
	printf("solver");
	char name[100];
	*result = sub_h.V;
	int error;
	GRBenv   *env   = NULL;
	GRBmodel *model = NULL;
	int qtd_labels_aux = sub_h.L - (sub_h.V - 1);
	int ind[qtd_labels_aux];
	double val[qtd_labels_aux];
	double val2[sub_h.L];
	int ind2[sub_h.L];
	int solcount;
	double *solution = (double*) malloc(qtd_labels_aux* sizeof(double));
  	double aux = sub_h.V;
	error = GRBloadenv(&env, "solver.log");
	if (error) errorr(env);
	error = GRBnewmodel(env, &model, "klsf", 0, NULL, NULL, NULL, NULL, NULL);
	if (error) errorr(env);
	

	for(int l = 0; l < sub_h.L; l++){
		if(l < qtd_labels_aux){
			ind[l] = l;
			val[l] = 1.0;
			solution[l] = 0.;
			sprintf(name, "x_%d", l);
			error = GRBaddvar(model, 0, NULL, NULL, 0.0 , 0.0, 1.0, GRB_BINARY, name);
		}else{
			sprintf(name, "xx_%d", l);	
			error = GRBaddvar(model, 0, NULL, NULL, 1.0 , 0.0, 1.0, GRB_BINARY, name);
		}
		val2[l] = sub_h.edges_per_label[l].qtd_edges;
		ind2[l] = l;
		if (error) errorr(env);
	}
	sprintf(name, "cores");
	error = GRBaddconstr(model, qtd_labels_aux, ind, val, GRB_EQUAL, sub_h.k, name);
	if (error) errorr(env);
  	sprintf(name, "aux");
  	error = GRBaddconstr(model, sub_h.L, ind2, val2, GRB_GREATER_EQUAL, sub_h.V - 1, name);
	if (error) errorr(env);
  	
  	error = GRBwrite(model, "model.rlp.gz");
  	if(error) errorr(env);	
	error = GRBsetcallbackfunc(model, calback_func, (void *) &sub_h);
	if (error) errorr(env);
	error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_LAZYCONSTRAINTS, 1);
	if (error) errorr(env);
	error = GRBsetintparam(GRBgetenv(model), "OutputFlag", 0);
	if (error) errorr(env);
	error = GRBsetdblparam(GRBgetenv(model), "TimeLimit", 4.0);
	if (error) errorr(env);
	error = GRBoptimize(model);
	if (error) errorr(env);
	error = GRBgetintattr(model, GRB_INT_ATTR_SOLCOUNT, &solcount);
  	if (error) errorr(env);
 	if(solcount){
		error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &aux);
		*result = aux;
  		if (error) errorr(env);
  		error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, qtd_labels_aux, solution);
  		if (error) errorr(env);
  	}else{
  		*result = V+1;
  	}


  	
	GRBfreemodel(model);
	GRBfreeenv(env);
  	return(solution);
}

void best_to_vns(VNS *vns, int k){
	for(int i = 0; i < k; i++){
		vns->solution_now[i] = vns->best_know[i];
	}
	for(int i =0; i < vns->qtd_labels_sub; i++){
		vns->mapped_label[vns->list_labels_sub[i]] = -1;
		vns->list_labels_sub[i] = -1;
	}
}

int main(int argc, char **argv){

	double result =0;
	double best = 0;


	srand(time(0));

	
	Problem G, sub_h;



	FILE *f_lista_arqvs;
	FILE *saida;
	f_lista_arqvs = fopen("lista_arquivos.txt", "r");
	
	
	if (f_lista_arqvs == NULL) goto QUIT;

	VNS vns;
	char s_nome_arquivo[100];
	int div10 = 0;
	int error = 0;
	int aux = 0;
	int it = 0;
	double *solution_solver = NULL;
	long  t_total_in_solver = 0;
	time_t t_tudo = time(NULL);
	time_t t_star_solver = time(NULL);
	int* sol_to_avaliation;
	while(fscanf(f_lista_arqvs, "%s\n", s_nome_arquivo) != EOF )
	{
		G = start_problem();
		div10 += 1;
		printf("%s\n", s_nome_arquivo);
		error = read_file(&G, s_nome_arquivo);
		if(error) goto QUIT;

		//generate_individuals(populacao, G);
		sol_to_avaliation = (int*) malloc(G.L * sizeof(int));
		for(int i = 0; i < G.L; i++)
			sol_to_avaliation[i] = 0;

		t_tudo = time(NULL);
		vns = start_alg(G.V, G.E, G.L, G.k);
		printf("%d", vns.qtd_fixed);
		it = 1;
		best = (double) G.V;
		result = best;
		while(time(NULL) - t_tudo < 200 && it < MAX_ITERATION){
			pertubation_and_fix_values(&vns, G.k, G.L);
			printa_infos(vns, G.V, G.E, G.L, G.k);
			sub_h = gera_sub_h(vns.mapped_label, vns.qtd_labels_sub, vns.qtd_fixed, 1, G);
			t_star_solver = time(NULL);
			printf("\n");
			solution_solver = generate_model_and_solve(sub_h, &result, G.V);
			printa_problem(sub_h);
			result = round(result);
			if(result < best && result > 0.1){
				from_solver_to_vns(&vns, solution_solver);
				best = result;
				for(int i = 0; i < G.L; i++)
					sol_to_avaliation[i] = 0;
				for(int i = 0; i < G.k; i ++)
					sol_to_avaliation[vns.best_know[i]] = 1;
				if(best == 1){
					break;
				}
				it = 1;
			}else{
				best_to_vns(&vns, G.k);
				it+=1;
			}
			printa_infos(vns, G.V, G.E, G.L, G.k);
			
			if(time(NULL) - t_star_solver > 2.5)
				vns.qtd_labels_sub -= G.k;
			
			
			t_total_in_solver += time(NULL) - t_star_solver;
			
			printf("\n%s\nbest: %f\nit: %d",s_nome_arquivo, best,it);
			finish_problem(&sub_h);
			vns.iteracao++;
		}

		saida = fopen("saida.txt", "a");
		fprintf(saida, "%s:%ld:%ld:%lf\n", s_nome_arquivo, time(NULL) - t_tudo, t_total_in_solver, best);
		fclose(saida);
		for(int i = 0; i < G.L; i++)
			sol_to_avaliation[i] = 0;
		for(int i = 0; i < G.k; i ++)
			sol_to_avaliation[vns.best_know[i]] = 1;
		printf("\n");
		aux = avalia_solution(sol_to_avaliation, G);
		printf("%f", result);
		if(aux != best){
			printf("\n%d %lf\n", aux, best);
			goto QUIT;
		}
		free(sol_to_avaliation);
		printf("\n%lf", result);
		finish_problem(&G); 
		finish_vns(&vns);
	}

QUIT:
	if(sol_to_avaliation != NULL)
		free(sol_to_avaliation);
	finish_problem(&G); 
	finish_vns(&vns);
	if (saida != NULL)
		fclose(saida);

	if (f_lista_arqvs != NULL)
		fclose(f_lista_arqvs);
}



