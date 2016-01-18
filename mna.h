#include "csparse.h"
#define TABLE_SIZE 32			//Hash table size
#define STRING_SIZE 32
#define EPS 1e-14

char **nodeIDs;
int g2_elements, id_counter;	//number of group 2 elements and max node id number
int nz_number; 			//Approximate number of non-zero elements
int nz_number_C;
int nz_number_G;

cs *trans_C;
cs *trans_G;
cs *comp_tr_C;
cs *comp_tr_G;
cs *BE_A;
cs *BE_temp_C;
cs *TR_A;
cs *TR_temp_GC;

css *S_trans;
csn *N_trans;

double **G_no_sparse;
double **C_no_sparse;
double **C_with_h;
bool x_flag; 

struct node_hash{
	char terminal[STRING_SIZE];
	int id;
	
	struct node_hash *next;
};

struct node_hash *hashtable[TABLE_SIZE];


//Node maping functions
int node_maping();
int node_maping_hash();
int hash_function(char *word);
void print_list(struct node_hash *head);

//MNA analysis functions
int mna_analysis();
int LU_factorization();
int cholesky_factorization();
int linearSolver();
int DC_sweep_solver();
int print_plotter();
int cgSolver();
int BiCGSolver();
void preconditioner(double *vectorZ, double *vectorR);
void matrixVectorMultiplication(double **matrix, double *vector, double *res);
double euclideanNorm(double *vector);
double innerProduct(double *vector1, double *vector2);

//MNA matrices and vectors
double **matrix_A;
double *vector_b;
double *vector_x;
double *vector_e;
double *vector_e_prev;
double *vector_temp;
double *vector_b_dupl;

//L diagonal vector
double *L_diag;

//Pivoting vector
int *P;	

void map_plot_nodes();

//Sparse analysis functions
int nz_approximation();
cs *create_triplet();
cs *compression_function(cs *A);
int sparse_LU(cs *C);
int sparse_cholesky(cs *C);
int sparse_CG(cs *C);
void sparse_preconditioner(cs *C, double *vectorZ, double *vectorR);
int sparse_BiCG(cs *C);


//Sparse transient analysis functions
void trans_init_sparse();
int nz_approximation_trans();
void create_trans_triplet();
void trans_compression_function();
void BE_sparseA(double h);
void sparse_tr_LU(cs *A);
void BE_sparse_b(double time, double final_time);
void trans_solve();
void TR_sparseA(double h);
void TR_sparse_b(double time, double final_time);
void trans_free_sparse();

//No sparse transient analysis functions
void trans_init_no_sparse();
void GC_calculation();
void BE_no_sparseA(double h);
void BE_no_sparse_b(double time, double final_time);
void TR_no_sparseA(double h);
void TR_no_sparse_b(double time, double final_time);
void trans_free_no_sparse();
