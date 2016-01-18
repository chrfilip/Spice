#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "mna.h"
#include "parser.h"

int initial_loc =0;
int p_vector_ind =0;

void map_plot_nodes() {
	
	struct plot_node *curr_node;
	struct element *curr;
	
	for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) {

		for(curr = root; curr != NULL; curr=curr->nxt_element) {
			
			if(strcmp(curr->terminal_a, curr_node->original_name) ==0) {
				curr_node->mapped_id = curr->terminalMNA_a;
				break;
			}
			else if(strcmp(curr->terminal_b, curr_node->original_name) ==0) {
				curr_node->mapped_id = curr->terminalMNA_b;
				break;
			}
		}
	}
}

int hash_function(char *word) {

	unsigned long hash = 5381;
	int c;

	while ((c = *word++))
		hash = ((hash << 5) + hash) + c;
	return hash % TABLE_SIZE;
}



int hash_table_add (char *terminal) {
	
	int position;
	struct node_hash *curr, *new_node;
	
	//Calculate position
	position = hash_function(terminal);
	
	//Check if name already exists
	for(curr = hashtable[position]; curr != NULL; curr=curr->next) {
		if(!strcmp(terminal, curr->terminal))
			return (curr->id);			
	}
	
	//Adds a new node name to the beginning of th list pointed to by head
	new_node = (struct node_hash *)malloc(sizeof( struct node_hash));
	if (new_node == NULL) {
		printf("Malloc Error! (hash table)\n");
		return 0;
	}
	
	strcpy(new_node->terminal, terminal);
	new_node->id = id_counter;
	new_node->next = hashtable[position]->next;
	hashtable[position]->next = new_node;
	id_counter++;
	printf("Node %s ---> key %d\n", terminal, id_counter);
	return(new_node->id);
}

int node_maping_hash() {
	
	struct element *curr;
	int i;
		
	g2_elements = 0;
	id_counter = 0;
		
	//Initialize hash table
	for (i=0; i<TABLE_SIZE; i++) {
		hashtable[i] = (struct node_hash *)malloc(sizeof(struct node_hash));
		hashtable[i]->next = NULL;
	}
	
	printf("Node maping using hash.\n\n");
	
	//Add new Elements
	for(curr = root; curr != NULL; curr=curr->nxt_element) {
		
		if( curr->element_type == 'V' || curr->element_type == 'L')
			g2_elements++;
		
		if( strcmp(curr->terminal_a, "0") != 0 )
			curr->terminalMNA_a = hash_table_add(curr->terminal_a);
		if( strcmp(curr->terminal_b, "0") != 0 )
			curr->terminalMNA_b = hash_table_add(curr->terminal_b);
		
	}

	//Test, prints netlist with maped node names
	printf("-------------------------------------\n");
	printf("The given netlist with hash node maping (-1 is ground)\n\n");
	for( curr=root; curr != NULL; curr = curr->nxt_element) {
		printf("%c%s %d %d %.2e\n", curr->element_type, curr->name, curr->terminalMNA_a, curr->terminalMNA_b, curr->value);
	}

	id_counter--;
	
	return 1;
}


void print_list(struct node_hash *head){
	struct node_hash *curr;
	
	curr = head->next;
	while (curr != NULL) {
		printf("%d ", curr->id);
		curr = curr->next;
		
	}
	printf("\n");
}


int mna_analysis() {

	struct element *curr;
	int matrix_dim = id_counter + 1 + g2_elements;
	int voltSources = 0, i=0;


	//Allocate memory for Matrix A, set to 0
	matrix_A = (double **)malloc(matrix_dim * sizeof(double*));
	if (matrix_A == NULL) {
		printf("MNA analysis: Matrix A allocation Error!\n");
		return 0;
	}
	for( i=0; i<matrix_dim; i++) {
		matrix_A[i] = (double *)calloc(matrix_dim, sizeof(double));
		if (matrix_A[i] == NULL) {
			printf("MNA analysis: Matrix A allocation Error!\n");
			return 0;
		}
	}
	
	//Allocate memory for vector b, set to 0
	vector_b = (double *)calloc( matrix_dim, sizeof(double));
	if (vector_b == NULL) {
		printf("MNA analysis: Vector b allocation Error!\n");
		return 0;
	}
	
	//Allocate memory for pivoting vector P
	P = (int *)malloc(matrix_dim * sizeof(int));
	if (P == NULL) {
		printf("MNA analysis: Vector P allocation Error!\n");
		return 0;
	}
	
	//Initialize vector P
	for (i=0; i<matrix_dim; i++)
		P[i] = i;
	
	//Calculate MNA system
	for( curr=root; curr != NULL; curr = curr->nxt_element) {
	  
	
		if (curr->element_type == 'R') {

			//Matrix A contribution
			if( curr->terminalMNA_a != -1)
				matrix_A[curr->terminalMNA_a][curr->terminalMNA_a] += 1.0 / (curr->value);
			if( curr->terminalMNA_b != -1)
				matrix_A[curr->terminalMNA_b][curr->terminalMNA_b] += 1.0 / (curr->value);

			if( (curr->terminalMNA_a != -1) && ( curr->terminalMNA_b != -1) ) {
				matrix_A[curr->terminalMNA_a][curr->terminalMNA_b] -= 1.0 / (curr->value);
				matrix_A[curr->terminalMNA_b][curr->terminalMNA_a] -= 1.0 / (curr->value);
			}

		} 
		else if (curr->element_type == 'I') {
			
			//Vector b contribution
			if (curr->terminalMNA_a != -1)			
				vector_b[curr->terminalMNA_a] -= curr->value;
			if (curr->terminalMNA_b != -1)			
				vector_b[curr->terminalMNA_b] += curr->value;
		}
		else if (curr->element_type == 'V') {

			//Matrix A contribution	
			if( curr->terminalMNA_a != -1) {
				matrix_A[curr->terminalMNA_a][id_counter + 1 + voltSources] += 1.0;
				matrix_A[id_counter + 1 + voltSources][curr->terminalMNA_a] += 1.0;
			}

			if( curr->terminalMNA_b != -1) {
				matrix_A[curr->terminalMNA_b][id_counter + 1 + voltSources] -= 1.0;
				matrix_A[id_counter + 1 + voltSources][curr->terminalMNA_b] -= 1.0;
			}

			//Vector b contribution
			vector_b[id_counter + 1 + voltSources] += curr->value;
			
			voltSources++;
		}
		else if (curr->element_type == 'L') {

			//Matrix A contribution	
			if( curr->terminalMNA_a != -1) {
				matrix_A[ curr->terminalMNA_a][id_counter + 1 + voltSources] += 1.0;
				matrix_A[id_counter + 1 + voltSources][curr->terminalMNA_a] += 1.0;
			}

			if( curr->terminalMNA_b != -1) {
				matrix_A[curr->terminalMNA_b][id_counter + 1 + voltSources] -= 1.0;
				matrix_A[id_counter + 1 + voltSources][curr->terminalMNA_b] -= 1.0;
			}
			
			vector_b[id_counter + 1 + voltSources] = 0;
			voltSources++;
		}
		
	}

	printf("Matrix A and vector b calculation complete!\n");
	
        x_flag = 1; // useful for no sparse transient analysis && DC_sweep

        if (trans_option) {
          vector_b_dupl = (double *) malloc(matrix_dim * sizeof(double)); 

          for (i = 0; i < matrix_dim; i++) {
            vector_b_dupl[i] = vector_b[i];
          }
        }      

	return 1;
}


int LU_factorization() {


        int matrix_dim = id_counter + 1 + g2_elements;

	int i, j, k, m, tmp_int;
	double x, *tmp_double;

	for(k=0; k<matrix_dim; k++){
		x = fabs(matrix_A[k][k]);
		m=k;
		// find the largest m for pivoting		
		for(i=k+1; i<matrix_dim; i++){
			if(fabs(matrix_A[i][k])>x) m=i;
			}
		if(matrix_A[m][k]==0){
			printf("Division by zero during LU factorization.\n");
			exit(1);
			}

		// make the pivot in matrix A and also in pivoting vector P			
		tmp_double = matrix_A[m];
		matrix_A[m] = matrix_A[k];
		matrix_A[k] = tmp_double;
		tmp_int = P[m];
		P[m] = P[k];
		P[k] = tmp_int;

		// make the LU calculations
		for(i=k+1; i<matrix_dim; i++){
			matrix_A[i][k] = matrix_A[i][k]/matrix_A[k][k]; 
			}
		for(i=k+1; i<matrix_dim; i++){
			for(j=k+1; j<matrix_dim; j++){
				matrix_A[i][j] = matrix_A[i][j] - matrix_A[i][k] * matrix_A[k][j];
				}
			}
		}

		if (x_flag) {

			//Vector x memory allocation
			vector_x = (double*)calloc(matrix_dim, sizeof(double));

                }
                x_flag = 0;       
		

	return 1;
	
}


int cholesky_factorization() {

	int i, k, j;
	int matrix_dim = id_counter + 1 + g2_elements;
	L_diag = (double*)malloc(matrix_dim * sizeof(double));
	double temp = 0.0;
	
	if(L_diag == NULL) {
		printf("L_diag allocation in cholesky error!\n");
		return 0;
	}
	
	for(k=0; k<matrix_dim; k++) {
		for(j=0; j<=k-1; j++)
			temp += pow(matrix_A[k][j], 2);
	
		temp = matrix_A[k][k] - temp;
		
		if(temp < 0) {
			printf("Error in cholesky function, matrix not SPD.\n");
			return 0;
		}
		
		L_diag[k] = matrix_A[k][k] = sqrt(temp); 

		temp = 0;
		for(i=k+1; i<matrix_dim; i++) {
			
			for(j=0; j<=k-1; j++)
				temp += matrix_A[i][j] * matrix_A[k][j];
			
			matrix_A[i][k] = (matrix_A[k][i] - temp) / L_diag[k];
			
			matrix_A[k][i] = matrix_A[i][k];		//filling symmetric elements

		}
		temp = 0.0;

	}
	
	return 1;
}

int linearSolver() {
	

        struct plot_node *curr_node; 
        int matrix_dim = id_counter + 1 + g2_elements;
	int k, j, i;
	double *y = (double*)calloc(matrix_dim,sizeof(double));

	// solve the L*y=b system
	y[0] = vector_b[P[0]];
	for(k=1; k<matrix_dim; k++){
		double sum;
		sum = 0.0;
		for(j=0; j<k; j++){
			sum += matrix_A[k][j]*y[j]; 
			}
		y[k] = vector_b[P[k]] - sum;
		}

	// solve the U*x=y system
	for(k=matrix_dim - 1; k>=0; k--){
		double sum;
		sum = 0.0;
		for(j=matrix_dim-1; j>k; j--){
			sum += matrix_A[k][j] * vector_x[j]; 
			}
		vector_x[k] = (y[k] - sum)/matrix_A[k][k];
		}

	
	
	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
			fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);	
	}
	else if(plot && trans_option){
		tran_output = fopen("tran_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
			fprintf(tran_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(tran_output);
	}
	else{
	  printf("x=\n");

	  for(i=0; i<matrix_dim; i++)
		printf("\t %f\n",vector_x[i]);
	}


	free(y);

	
	return 1;	

}


int DC_sweep_solver() {
	
	int matrix_dim = id_counter + 1 + g2_elements;
	struct element *curr; 
	char tempstr[STRING_SIZE] = "", temp[2];
	int sourcectr = 0;
	int i=0;
	double j;
	
	for(curr = root; curr != NULL; curr=curr->nxt_element) {
			
		temp[0] = curr->element_type;
		if(temp[0] == 'V')
			sourcectr++;
		temp[1] = '\0';
		strcpy(tempstr, temp );
		strcat(tempstr, curr->name);
		if(strcmp(tempstr, dc_source) == 0) {
			break;
		}
		i++;
	}

	dc_output = fopen("dc_output.txt", "w");
	fclose(dc_output);
	for(j=dc_start, i=0; j<=dc_end+0.001; j+=dc_step, i++) {
		
		vector_b[i] = j;
		dc_output = fopen("dc_output.txt", "a");
		fprintf(dc_output, "Iteration: %d\n", i+1);
		fclose(dc_output);
		if(sparse_option){
		  if(SPD_option && ITER_option){
		    sparse_CG(compression_function(create_triplet()));
		  }
		  else if(SPD_option && !ITER_option){
		    sparse_cholesky(compression_function(create_triplet()));
		  }
		  else if(!SPD_option && ITER_option){
		    sparse_BiCG(compression_function(create_triplet()));
		  }
		  else if(!SPD_option && !ITER_option){
		    sparse_LU(compression_function(create_triplet()));
		  }
		}
		else{
		  if(SPD_option && ITER_option){
		    cgSolver();
		  }
		  else if(SPD_option && !ITER_option){
		    linearSolver();
		  }
		  else if(!SPD_option && ITER_option){
		    BiCGSolver();
		  }
		  else if(!SPD_option && !ITER_option){
		  linearSolver();
		  }
		}
		
	}
	if(plot && DC_sweep){
	  
	  plot = DC_sweep = 0;
	  
	  if(sparse_option){
		  if(SPD_option && ITER_option){
		    sparse_CG(compression_function(create_triplet()));
		  }
		  else if(SPD_option && !ITER_option){
		    sparse_cholesky(compression_function(create_triplet()));
		  }
		  else if(!SPD_option && ITER_option){
		    sparse_BiCG(compression_function(create_triplet()));
		  }
		  else if(!SPD_option && !ITER_option){
		    sparse_LU(compression_function(create_triplet()));
		  }
		}
		else{
		  if(SPD_option && ITER_option){
		    cgSolver();
		  }
		  else if(SPD_option && !ITER_option){
		    linearSolver();
		  }
		  else if(!SPD_option && ITER_option){
		    BiCGSolver();
		  }
		  else if(!SPD_option && !ITER_option){
		  linearSolver();
		  }
		}
	  
	}

	return 1;
}


int cgSolver() {

	struct plot_node *curr_node;
	int k;
	int i,iter = 0;
	int matrix_dim = id_counter + 1 + g2_elements;
	double temp, rho, beta, rho1, alpha;
	double *vector_r = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_z = (double*)calloc(matrix_dim, sizeof(double));
	double *matrixVectorMul = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_p = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_q = (double*)calloc(matrix_dim, sizeof(double));

	vector_x = (double*)calloc(matrix_dim, sizeof(double));
	
	//Calculate Ax
	matrixVectorMultiplication(matrix_A, vector_x, matrixVectorMul);
	
	//Calculate r = b - Ax
	for(i=0; i< matrix_dim; i++) 
		vector_r[i] = vector_b[i] - matrixVectorMul[i];
	
	temp = euclideanNorm(vector_b);
	if(temp == 0.0)
		temp = 1.0;
	
	while( euclideanNorm(vector_r) / temp > itol && iter < matrix_dim) {

		iter++;

		preconditioner(vector_z, vector_r);
		
		rho = innerProduct(vector_r, vector_z);
		
		if(iter == 1) {
			for(i=0; i<matrix_dim; i++)
				vector_p[i] = vector_z[i];
		} 
		else {
			beta = rho / rho1;
			for(i=0; i<matrix_dim; i++)
				vector_p[i] = vector_z[i] + beta * vector_p[i];
		}
		
		rho1 = rho;
		
		matrixVectorMultiplication(matrix_A, vector_p, vector_q);
		
		alpha = rho / innerProduct(vector_p, vector_q);
		
		for(i=0; i<matrix_dim; i++) {
				vector_x[i] = vector_x[i] + alpha * vector_p[i];
				vector_r[i] = vector_r[i] - alpha * vector_q[i];
		}

	}
	
	
	
	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);	
	}
	else{
	  printf("x=\n");
	
	  for(i=0; i<matrix_dim; i++)
		printf("\t %f\n",vector_x[i]);

	}
	
	free(vector_p);
	free(vector_r);
	free(vector_z);
	free(vector_q);
	free(matrixVectorMul);
	
	free(vector_x);
	
	return 1;
}

int BiCGSolver() {

	struct plot_node *curr_node;
	int k;
	int i,iter = 0, j;
	int matrix_dim = id_counter + 1 + g2_elements;
	double temp, rho, beta, rho1, alpha, omega;
	double *vector_r = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_rtld = (double*)calloc(matrix_dim, sizeof(double));

	double *vector_z = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_ztld = (double*)calloc(matrix_dim, sizeof(double));
	
	double *matrixVectorMul = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_p = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_ptld = (double*)calloc(matrix_dim, sizeof(double));
	
	double *vector_q = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_qtld = (double*)calloc(matrix_dim, sizeof(double));
	
	double **transposedA = (double **)malloc(matrix_dim * sizeof(double*));
	if (transposedA == NULL) {
		printf("Matrix transposed A allocation Error!\n");
		return 0;
	}
	for(i=0; i<matrix_dim; i++) {
		transposedA[i] = (double *)calloc(matrix_dim, sizeof(double));
	}

	vector_x = (double*)calloc(matrix_dim, sizeof(double));
	
	//Calculate Ax
	matrixVectorMultiplication(matrix_A, vector_x, matrixVectorMul);
	
	//Calculate r = r_tld = b - Ax
	for(i=0; i< matrix_dim; i++) 
		vector_r[i] = vector_rtld[i] = vector_b[i] - matrixVectorMul[i];
	
	temp = euclideanNorm(vector_b);
	if(temp == 0.0)
		temp = 1.0;
	
	while( euclideanNorm(vector_r) / temp > itol && iter < matrix_dim) {
		
		iter++;
		
		preconditioner(vector_z, vector_r);
		preconditioner(vector_ztld, vector_rtld);
		
		rho = innerProduct(vector_z, vector_rtld);
		
		if( fabs(rho) < EPS) {
			return 0;
		}
	
		if(iter == 1) {
			
			for(i=0; i<matrix_dim; i++) {
				vector_p[i] = vector_z[i];
				vector_ptld[i] = vector_ztld[i];
			}
		} 
		else {
			
			beta = rho / rho1;
			
			for(i=0; i<matrix_dim; i++) {
				vector_p[i] = vector_z[i] + beta * vector_p[i];
				vector_ptld[i] = vector_ztld[i] + beta * vector_ptld[i];
			}
		}
		
		rho1 = rho;
		
		matrixVectorMultiplication(matrix_A, vector_p, vector_q);

		//Calculate transposed A
		for(i=0; i<matrix_dim; i++) {
			for(j=0; j<matrix_dim; j++)
				transposedA[i][j] = matrix_A[j][i];
		}
		
		matrixVectorMultiplication(transposedA, vector_ptld, vector_qtld);

		omega = innerProduct(vector_ptld, vector_q);
		
		if( fabs(omega)<EPS){
			return 0;
		}
		alpha = rho / omega;
		
		for(i=0; i<matrix_dim; i++) {
			vector_x[i] = vector_x[i] + alpha * vector_p[i];
			vector_r[i] = vector_r[i] - alpha * vector_q[i];
			vector_rtld[i] = vector_rtld[i] - alpha * vector_qtld[i];

		}

	}
	
	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);	
	}
	else{
	  printf("x=\n");
	
	  for(i=0; i<matrix_dim; i++)
		printf("\t %f\n",vector_x[i]);

	}
	
	free(vector_p);
	free(vector_r);
	free(vector_z);
	free(vector_q);
	free(vector_ptld);
	free(vector_rtld);
	free(vector_ztld);
	free(vector_qtld);

        for(i=0; i<matrix_dim; i++) { 
	  free(transposedA[i]);
        } 

        free(transposedA); 

	free(matrixVectorMul);

        free(vector_x);
	
	return 1;

}

double euclideanNorm(double *vector) {
	
	int i;
	double temp = 0.0;
	int matrix_dim = id_counter + 1 + g2_elements;
	
	for(i=0; i<matrix_dim; i++)
		temp += pow(vector[i],2);
	
	return sqrt(temp);
}

void matrixVectorMultiplication(double **matrix, double *vector, double *res) {

	int i,j;
	double temp = 0.0;
	int matrix_dim = id_counter + 1 + g2_elements;
	
	for (i=0; i<matrix_dim; i++) {
		for(j=0; j<matrix_dim; j++)
			temp += matrix[i][j] * vector[j];
		res[i] = temp;
		temp = 0.0;
	}
	
	return;
}

double innerProduct( double *vector1, double *vector2) {
	
	int i;
	int matrix_dim = id_counter + 1 + g2_elements;
	double result = 0.0;
	
	for(i=0; i<matrix_dim; i++)
		result += vector1[i] * vector2[i];
	
	return result;
	
}

void preconditioner(double *vectorZ, double *vectorR) {
	
	int i;
	int matrix_dim = id_counter + 1 + g2_elements;
	double temp;
	
	for(i=0; i<matrix_dim; i++) {
		
		temp = matrix_A[i][i];
		if( temp == 0.0)
			temp = 1.0;
		
		vectorZ[i] =( 1.0/temp ) * vectorR[i];
		
	}
	
	return;
}

void sparse_preconditioner(cs *C, double *vectorZ, double *vectorR) {
	
	int j = 0, k;
	int matrix_dim = id_counter + 1 + g2_elements;

        int diag_elem = 0;
        int check = 0; 
   
        while ( C->p[j] < C->p[matrix_dim] ) {          

          for (k = C->p[j]; k < C->p[j+1]; k++) {

            if ( C->i[k] == diag_elem ) {
              vectorZ[j] =( 1.0/C->x[k] ) * vectorR[j];
              diag_elem++; 
              check = 1;
              break;
            }
             
          }

          if ( check != 1 ) {
            vectorZ[j] = 0;
            diag_elem++;
          }   
          
	  check = 0;      
          j++;
        }          
         	
	return;
	
}

int nz_approximation() {

	struct element *curr;
	int voltSources = 0;
	int matrix_dim = id_counter + 1 + g2_elements;
	int i;
	
    nz_number = 0; 

	//Allocate memory for vector b, set to 0
	vector_b = (double *)calloc( matrix_dim, sizeof(double) );
	if (vector_b == NULL) {
		printf("Sparse analysis: Vector b allocation Error!\n");
		return 0;
	}
	
	for (curr = root; curr != NULL; curr = curr->nxt_element) {
		
		if (curr->element_type == 'R') {

			//Resistance contribution 
			if( curr->terminalMNA_a != -1)
				nz_number += 1;
			if( curr->terminalMNA_b != -1)
				nz_number += 1;
			if( (curr->terminalMNA_a != -1) && ( curr->terminalMNA_b != -1) ) 
				nz_number += 2;

		} 
		else if (curr->element_type == 'V') {

			//Voltage source contribution 
			if( curr->terminalMNA_a != -1) 
				nz_number += 2;
			if( curr->terminalMNA_b != -1) 
				nz_number += 2;

			//Vector b contribution
			vector_b[id_counter + 1 + voltSources] += curr->value;

                        voltSources++;
		}
		else if (curr->element_type == 'I') {
			
			//Vector b contribution
			if (curr->terminalMNA_a != -1)			
				vector_b[curr->terminalMNA_a] -= curr->value;
			if (curr->terminalMNA_b != -1)			
				vector_b[curr->terminalMNA_b] += curr->value;
		}
		else if (curr->element_type == 'L') {

			//L contribution 
			if( curr->terminalMNA_a != -1) 
				nz_number += 2;
			if( curr->terminalMNA_b != -1) 
				nz_number += 2;		

                        voltSources++;
		}
		
	}

        if (trans_option) {
          vector_b_dupl = (double *) malloc(matrix_dim * sizeof(double)); 

          for (i = 0; i < matrix_dim; i++) {
            vector_b_dupl[i] = vector_b[i];
          }
        }       

	printf("Non-zero calculation completed. %d\n", nz_number);

	return 1;
}

cs  *create_triplet() {

	int matrix_dim = id_counter + 1 + g2_elements;
	int voltageSources = 0;
	int k = 0;
	struct element *curr;
	
	cs *A = cs_spalloc(matrix_dim,matrix_dim,nz_number,1,1);
	
	for(curr=root; curr != NULL; curr = curr->nxt_element) {
		
		if (curr->element_type == 'R') {

			//Resistance contribution 
			if( curr->terminalMNA_a != -1) {
				A->i[k] = curr->terminalMNA_a;
				A->p[k] = curr->terminalMNA_a;
				A->x[k] = 1.0/curr->value;
				k++;
			}

			if( curr->terminalMNA_b != -1) {
				A->i[k] = curr->terminalMNA_b;
				A->p[k] = curr->terminalMNA_b;
				A->x[k] = 1.0/curr->value;
				k++;
                        }

			if( ( curr->terminalMNA_a != -1) && ( curr->terminalMNA_b != -1) ) {
				A->i[k] = curr->terminalMNA_a;
				A->p[k] = curr->terminalMNA_b;
				A->x[k] = -1.0/curr->value;
				k++;

				A->i[k] = curr->terminalMNA_b;
				A->p[k] = curr->terminalMNA_a;
				A->x[k] = -1.0/curr->value;
				k++;
			}

		} 
		else if (curr->element_type == 'V') {

			//Voltage source contribution 
			if( curr->terminalMNA_a != -1) {				
				A->i[k] = curr->terminalMNA_a;
				A->p[k] = id_counter + 1 + voltageSources;
				A->x[k] = 1.0;
				k++;

				A->i[k] = id_counter + 1 + voltageSources;
				A->p[k] = curr->terminalMNA_a;
				A->x[k] = 1.0;
				k++;
			}

			if( curr->terminalMNA_b != -1) {
				A->i[k] = curr->terminalMNA_b;
				A->p[k] = id_counter + 1 + voltageSources;
				A->x[k] = -1.0;
				k++;

				A->i[k] = id_counter + 1 + voltageSources;
				A->p[k] = curr->terminalMNA_b;
				A->x[k] = -1.0;
				k++;
			}

			voltageSources++;
		}
		else if (curr->element_type == 'L') {

			//L contribution 
			if( curr->terminalMNA_a != -1) {
				A->i[k] = curr->terminalMNA_a;
				A->p[k] = id_counter + 1 + voltageSources;
				A->x[k] = 1.0;
				k++;

				A->i[k] = id_counter + 1 + voltageSources;
				A->p[k] = curr->terminalMNA_a;
				A->x[k] = 1.0;
				k++;
			}

			if( curr->terminalMNA_b != -1) {
				A->i[k] = curr->terminalMNA_b;
				A->p[k] = id_counter + 1 + voltageSources;
				A->x[k] = -1.0;
				k++;

				A->i[k] = id_counter + 1 + voltageSources;
				A->p[k] = curr->terminalMNA_b;
				A->x[k] = -1.0;
				k++;	
			}

			voltageSources++;
		}
		
	}

	A->nz = k;

	printf("Triplet creation process completed.\n");

	return A;
}

cs *compression_function(cs *A) {
	
	int res;

	cs *C = cs_compress(A);

	cs_spfree(A);
	
	res = cs_dupl(C);
	
	if (!res) {
		printf("cs_dupl unsuccessful\n");
		return NULL;
	}

	return C; 	
}

int sparse_LU(cs *C) {
	
	css *S;
	csn *N;
	int matrix_dim = id_counter + 1 + g2_elements;
	int k;  
	struct plot_node *curr_node;
	
	//Vector x memory allocation
	vector_x = (double *)calloc( matrix_dim, sizeof(double) );
	if (vector_x == NULL) {
		printf("Sparse analysis: Vector x allocation Error!\n");
		return 0;
	}
	
	S = cs_sqr(2,C,0);
	N = cs_lu(C,S,1);
	cs_spfree(C);
	
	cs_ipvec(N->pinv,vector_b,vector_x,matrix_dim);
	cs_lsolve(N->L,vector_x);
	cs_usolve(N->U,vector_x);
	cs_ipvec(S->q,vector_x,vector_b,matrix_dim);
	
	cs_sfree(S);
	cs_nfree(N);
	
	
	
	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);
	}
	else{
	  printf("x=\n");
	    for(k=0; k<matrix_dim; k++)
	      printf("\t %f\n",vector_x[k]);
	 
	}
	
	free(vector_x);
	
	return 1;
}

int sparse_cholesky(cs *C) {

	css *S;
	csn *N;
	int matrix_dim = id_counter + 1 + g2_elements;
        int k;  
	struct plot_node *curr_node;

	//Vector x memory allocation
	vector_x = (double *)calloc( matrix_dim, sizeof(double) );
	if (vector_x == NULL) {
		printf("Sparse analysis: Vector x allocation Error!\n");
		return 0;
	}

        S = cs_schol(1,C);
        N = cs_chol(C,S);
        cs_spfree(C);

        cs_ipvec(S->pinv,vector_b,vector_x,matrix_dim);
        cs_lsolve(N->L,vector_x);
        cs_ltsolve(N->L,vector_x);
        cs_pvec(S->pinv,vector_x,vector_b,matrix_dim); 

        cs_sfree(S);
        cs_nfree(N);

	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);
	}
	else{
	  printf("x=\n");
	    for(k=0; k<matrix_dim; k++)
	      printf("\t %f\n",vector_x[k]);
	 
	}
        free(vector_x);

	return 1;
}

int sparse_CG(cs *C) {

	struct plot_node *curr_node;
	int k;
	int i,iter = 0;
	int matrix_dim = id_counter + 1 + g2_elements;
	double temp, rho, beta, rho1, alpha;
	double *vector_r = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_z = (double*)calloc(matrix_dim, sizeof(double));
	double *matrixVectorMul = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_p = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_q = (double*)calloc(matrix_dim, sizeof(double));

	vector_x = (double*)calloc(matrix_dim, sizeof(double));
	
	//Calculate Cx
	cs_gaxpy(C, vector_x, matrixVectorMul);
	
	//Calculate r = b - Cx
	for(i=0; i< matrix_dim; i++) 
		vector_r[i] = vector_b[i] - matrixVectorMul[i];
	
	temp = euclideanNorm(vector_b);
	if(temp == 0.0)
		temp = 1.0;
	
	while( euclideanNorm(vector_r) / temp > itol && iter < matrix_dim) {

		iter++;

		sparse_preconditioner(C, vector_z, vector_r); 
		
		rho = innerProduct(vector_r, vector_z);
		
		if(iter == 1) {
			for(i=0; i<matrix_dim; i++)
				vector_p[i] = vector_z[i];
		} 
		else {
			beta = rho / rho1;
			for(i=0; i<matrix_dim; i++)
				vector_p[i] = vector_z[i] + beta * vector_p[i];
		}
		
		rho1 = rho;
		
                for (i = 0; i < matrix_dim; i++) {
                  vector_q[i] = 0.0;  
                }
		cs_gaxpy(C, vector_p, vector_q);
		
		alpha = rho / innerProduct(vector_p, vector_q);
		
		for(i=0; i<matrix_dim; i++) {
				vector_x[i] = vector_x[i] + alpha * vector_p[i];
				vector_r[i] = vector_r[i] - alpha * vector_q[i];
		}

	}
	
	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);
	}
	else{
	  printf("x=\n");
	    for(k=0; k<matrix_dim; k++)
	      printf("\t %f\n",vector_x[k]);
	 
	}
	
	free(vector_p);
	free(vector_r);
	free(vector_z);
	free(vector_q);
	free(matrixVectorMul);
        cs_spfree(C);

        free(vector_x);
	
	return 1;

}


int sparse_BiCG(cs *C) {

  	struct plot_node *curr_node;
	int k;
	int i,iter = 0, j;
	int matrix_dim = id_counter + 1 + g2_elements;
	double temp, rho, beta, rho1, alpha, omega;
	double *vector_r = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_rtld = (double*)calloc(matrix_dim, sizeof(double));

	double *vector_z = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_ztld = (double*)calloc(matrix_dim, sizeof(double));
	
	double *matrixVectorMul = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_p = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_ptld = (double*)calloc(matrix_dim, sizeof(double));
	
	double *vector_q = (double*)calloc(matrix_dim, sizeof(double));
	double *vector_qtld = (double*)calloc(matrix_dim, sizeof(double));


	vector_x = (double*)calloc(matrix_dim, sizeof(double));
	
	//Calculate Cx
	cs_gaxpy(C, vector_x, matrixVectorMul);
	
	//Calculate r = r_tld = b - Ax
	for(i=0; i< matrix_dim; i++) 
		vector_r[i] = vector_rtld[i] = vector_b[i] - matrixVectorMul[i];
	
	temp = euclideanNorm(vector_b);
	if(temp == 0.0)
		temp = 1.0;
	
	while( euclideanNorm(vector_r) / temp > itol && iter < matrix_dim) {
		
		iter++;
		
		sparse_preconditioner(C, vector_z, vector_r);
		sparse_preconditioner(C, vector_ztld, vector_rtld);
		
		rho = innerProduct(vector_z, vector_rtld);
		
		if( fabs(rho) < EPS) {
			return 0;
		}
	
		if(iter == 1) {
			
			for(i=0; i<matrix_dim; i++) {
				vector_p[i] = vector_z[i];
				vector_ptld[i] = vector_ztld[i];
			}
		} 
		else {
			
			beta = rho / rho1;
			
			for(i=0; i<matrix_dim; i++) {
				vector_p[i] = vector_z[i] + beta * vector_p[i];
				vector_ptld[i] = vector_ztld[i] + beta * vector_ptld[i];
			}
		}
		
		rho1 = rho;
		for (i = 0; i < matrix_dim; i++) {
                  vector_q[i] = 0.0;  
                }
		cs_gaxpy(C, vector_p, vector_q);


                for (i = 0; i < matrix_dim; i++) {
                  for (j = C->p[i]; j < C->p[i+1]; j++)
                    vector_qtld[i] += C->x[j] * vector_ptld[C->i[j]];
                }    

		omega = innerProduct(vector_ptld, vector_q);
		
		if( fabs(omega)<EPS){
			return 0;
		}
		alpha = rho / omega;
		
		for(i=0; i<matrix_dim; i++) {
			vector_x[i] = vector_x[i] + alpha * vector_p[i];
			vector_r[i] = vector_r[i] - alpha * vector_q[i];
			vector_rtld[i] = vector_rtld[i] - alpha * vector_qtld[i];

		}

	}

	if(plot && DC_sweep) {
		dc_output = fopen("dc_output.txt", "a");
		for(curr_node = plot_root; curr_node != NULL; curr_node = curr_node->next) 
		  fprintf(dc_output, "Node %s: %lf\n", curr_node->original_name,vector_x[curr_node->mapped_id]);
		fclose(dc_output);
	}
	else{
	  printf("x=\n");
	    for(k=0; k<matrix_dim; k++)
	      printf("\t %f\n",vector_x[k]);
	 
	}
	
	free(vector_p);
	free(vector_r);
	free(vector_z);
	free(vector_q);
	free(vector_ptld);
	free(vector_rtld);
	free(vector_ztld);
	free(vector_qtld);
	free(matrixVectorMul);
        cs_spfree(C);

        free(vector_x); 
	
	return 1;

}





 //Given 2 points (x0,y0) and (x1,y0) and an x coordinate of a third point, it calculates the y coordinate of this
 //third point using Linear Interpolation. Returns this y coordinate.
double linearInterpolation(double x0, double x1, double y0, double y1, double x) {
	return ((y1-y0)/(x1-x0)) * (x-x0) + y0;
}



 //Computes the [transient_spec] EXP(i1 i2 td1 tc1 td2 tc2). The parameters time and final_time are required for
 //finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries,
 //the function returns a negative value.

double EXP_func(double i1, double i2, double td1, double tc1, double td2, double tc2, double time, double final_time) {

	double pow1, pow2;

	if (time>=0.0 && time<=td1) return i1;
	else if (time>td1 && time<=td2) {
		pow1 = (-1.0)*(time-td1)/tc1;
		return (i1+(i2-i1)*(1.0-exp(pow1)));
	}
	else if (time>td2 && time<=final_time) {
		pow1 = (-1.0)*(time-td1)/tc1;
		pow2 = (-1.0)*(time-td2)/tc2;
		return (i1+(i2-i1)*(exp(pow2)-exp(pow1)));
	}    
	else{
		printf("invalid time in SIN calculation\n");
		exit(1);
	}

}



//This function computes the [transient_spec] SIN(i1 ia fr td df ph). The parameters time and final_time are required
//for finding which expression to be used. Returns the computed value. If the parameter time is out of time boundaries
//the function returns a negative value.

double SIN_func(double i1, double ia, double fr, double td, double df, double ph, double time, double final_time) {

	double pi = 3.1459;
	double pow;

	if (time>=0.0 && time<=td) return (i1+ia*sin(2*pi*ph/360));
	else if (time>td && time<=final_time) { 
		pow = (-1.0) * (time-td) * df;  
		return (i1+ia*sin(2*pi*fr*(time-td)+2*pi*ph/360)*exp(pow));
	}
	else {
 		printf("invalid time in EXP calculation\n");
		exit(1);
	}
}



//This function computes the [transient_spec] PULSE(i1 i2 td tr tf pw per). The parameters time and final_time are
//required for finding which expression to be used. It returns the computed value. If the parameter time is out of
//time boundaries, the function returns a negative value.

double PULSE_func(double i1, double i2, double td, double tr, double tf, double pw, double per, double time, double final_time) {

	int k;

	k = (int)(time/per);
	if (time>=k*per && time<(td+k*per)) return i1;
	else if (time>=(td+k*per) && time<(td+tr+k*per)) return linearInterpolation(td+k*per, td+tr+k*per, i1, i2, time);
	else if (time>=(td+tr+k*per) && time<(td+tr+pw+k*per)) return i2;
	else if (time>=(td+tr+pw+k*per) && time<(td+tr+pw+tf+k*per)) return linearInterpolation(td+tr+pw+k*per, td+tr+pw+tf+k*per, i2, i1, time);
	else if (time>=(td+tr+pw+tf+k*per) && time<=(td+per+k*per)) return i1;
	else printf("invalid time in PULSE calculation\n");

	exit(1);	

}



//This function computes the [transient_spec] PWL(t1 v1 t2 v2 ... tn vn). The parameters time and final_time are
//required for finding which expression to be used. It returns the computed value. If the parameter time is out of
//time boundaries, the function returns a negative value.

double PWL_func(int n, double *t, double *i, double time, double final_time) {

	int k;

	for(k=0; k<n; k++)
		if(t[k]<=time && time<=t[k+1])   
			return linearInterpolation(t[k], t[k+1], i[k], i[k+1], time);
	printf("invalid time in PWL calculation\n");
	exit(1);		

}



void trans_free_sparse() {
	free(vector_e);
        free(vector_e_prev);
    	free(vector_temp);

}

void trans_init_no_sparse() {
	
	int i;
	int matrix_dim = id_counter + 1 + g2_elements;
	
	// Allocate memory for Matrix G, set to 0
	G_no_sparse = (double **)malloc(matrix_dim * sizeof(double*));
	if (G_no_sparse == NULL) {
		printf("Transient analysis: Matrix G allocation Error!\n");
		exit(1);
	}
	for (i = 0; i < matrix_dim; i++) {
		G_no_sparse[i] = (double *)calloc(matrix_dim, sizeof(double));
		if (G_no_sparse[i] == NULL) {
			printf("Transient analysis: Matrix G allocation Error!\n");
			exit(1);
		}
	}
	
	// Allocate memory for Matrix C, set to 0
	C_no_sparse = (double **)malloc(matrix_dim * sizeof(double*));
	if (C_no_sparse == NULL) {
		printf("Transient analysis: Matrix C allocation Error!\n");
		exit(1);
	}
	for (i = 0; i < matrix_dim; i++) {
		C_no_sparse[i] = (double *)calloc(matrix_dim, sizeof(double));
		if (C_no_sparse[i] == NULL) {
			printf("Transient analysis: Matrix C allocation Error!\n");
			exit(1);
		}
	}
	
	// Allocate memory for Matrix C_with_h, set to 0
	C_with_h = (double **)malloc(matrix_dim * sizeof(double*));
	if (C_with_h == NULL) {
		printf("Transient analysis: Matrix C_with_h allocation Error!\n");
		exit(1);
	}
	for (i = 0; i < matrix_dim; i++) {
		C_with_h[i] = (double *)calloc(matrix_dim, sizeof(double));
		if (C_with_h[i] == NULL) {
			printf("Transient analysis: Matrix C_with_h allocation Error!\n");
			exit(1);
		}
	}
	
	// Initialize vector P appropriately
	for (i = 0; i < matrix_dim; i++) { 
		P[i] = i;
	} 
	
	x_flag = 0;
	
	vector_e = (double *) malloc(matrix_dim * sizeof(double));
	vector_e_prev = (double *) malloc(matrix_dim * sizeof(double));  
	
}


void GC_calculation() {

	struct element *curr;

	int voltageSources = 0;
	

	for (curr = root; curr != NULL; curr = curr->nxt_element) {
		
		if (curr->element_type == 'R') {
			
			// Resistance contribution 
			if (curr->terminalMNA_a != -1) {
				G_no_sparse[curr->terminalMNA_a][curr->terminalMNA_a] += 1.0/curr->value;
			}
			
			if (curr->terminalMNA_b != -1) {
				G_no_sparse[curr->terminalMNA_b][curr->terminalMNA_b] += 1.0/curr->value;	
			}
			
			if ( (curr->terminalMNA_a != -1) && (curr->terminalMNA_b != -1) ) {
				G_no_sparse[curr->terminalMNA_a][curr->terminalMNA_b] -= 1.0/curr->value;	
				G_no_sparse[curr->terminalMNA_b][curr->terminalMNA_a] -= 1.0/curr->value; 
			}
			
		} 
		else if (curr->element_type == 'V') {
			
			// Voltage source contribution 
			if (curr->terminalMNA_a != -1) {
				G_no_sparse[curr->terminalMNA_a][id_counter + 1 + voltageSources] += 1.0;					
				G_no_sparse[id_counter + 1 + voltageSources][curr->terminalMNA_a] += 1.0;
			}
			
			if (curr->terminalMNA_b != -1) {
				G_no_sparse[curr->terminalMNA_b][id_counter + 1 + voltageSources] -= 1.0;
				G_no_sparse[id_counter + 1 + voltageSources][curr->terminalMNA_b] -= 1.0;
			}
			
			voltageSources++;
			
		}
		else if (curr->element_type == 'L') {
			
			// L contribution 
			if (curr->terminalMNA_a != -1) {
				G_no_sparse[curr->terminalMNA_a][id_counter + 1 + voltageSources] += 1.0; 
				G_no_sparse[id_counter + 1 + voltageSources][curr->terminalMNA_a] += 1.0; 
			}
			
			if (curr->terminalMNA_b != -1) {
				G_no_sparse[curr->terminalMNA_b][id_counter + 1 + voltageSources] -= 1.0;
				G_no_sparse[id_counter + 1 + voltageSources][curr->terminalMNA_b] -= 1.0; 
			}
			
			if ( (curr->terminalMNA_a != -1) && (curr->terminalMNA_b != -1) ) {
				C_no_sparse[id_counter + 1 + voltageSources][id_counter + 1 + voltageSources] -= curr->value;
			}
			
			voltageSources++;         
			
		}
		
	}
	
}


//This function calculates the matrix A of the transient no sparse linear system, using Backward Euler method

void BE_no_sparseA(double h) {

	int i, j;
        int matrix_dim = id_counter + 1 + g2_elements;

        for (i = 0; i < matrix_dim; i++) {
          for (j = 0; j < matrix_dim; j++) {  
            C_with_h[i][j] = C_no_sparse[i][j] * 1.0/h;
            matrix_A[i][j] = G_no_sparse[i][j] + C_with_h[i][j];  
          }
        }
  
}



//This function calculates vector b of the no sparse transient linear system, using Backward Euler method

void BE_no_sparse_b(double time, double final_time) {

	int matrix_dim = id_counter + 1 + g2_elements;
        int voltageSources = 0;
 
	struct element *curr;

	int i;
	float value;   

        for (i = 0; i < matrix_dim; i++) {
          vector_e[i] = 0.0; 
        }
	
        // calculate e(t) which changes at every step of transient analysis  
	for (curr = root; curr != NULL; curr = curr->nxt_element) {

		if (curr->trans_src != NULL) {

			if (curr->trans_src->type == 1) {
				value = EXP_func(curr->trans_src->trans_structs_union.trans_EXP.i1, curr->trans_src->trans_structs_union.trans_EXP.i2,
                                                 curr->trans_src->trans_structs_union.trans_EXP.td1, curr->trans_src->trans_structs_union.trans_EXP.tc1,
                                                 curr->trans_src->trans_structs_union.trans_EXP.td2, curr->trans_src->trans_structs_union.trans_EXP.tc2, 
                                                 time, final_time); 
                        }              
 			else if (curr->trans_src->type == 2) { 
				value = SIN_func(curr->trans_src->trans_structs_union.trans_SIN.i1, curr->trans_src->trans_structs_union.trans_SIN.ia,
                                                 curr->trans_src->trans_structs_union.trans_SIN.fr, curr->trans_src->trans_structs_union.trans_SIN.td,
                                                 curr->trans_src->trans_structs_union.trans_SIN.df, curr->trans_src->trans_structs_union.trans_SIN.ph,
                                                 time, final_time);
                        }
                        else if (curr->trans_src->type == 3) { 
				value = PULSE_func(curr->trans_src->trans_structs_union.trans_PULSE.i1, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.i2, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.td, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.tr, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.tf, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.pw, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.per, 
                                                   time, final_time);

                        }
                        else if (curr->trans_src->type == 4) {
 				value = PWL_func(curr->trans_src->trans_structs_union.trans_PWL.n, 
                                                 curr->trans_src->trans_structs_union.trans_PWL.t, 
                                                 curr->trans_src->trans_structs_union.trans_PWL.i, 
                                                 time, final_time);
                        }
                        
		}   
		else {
                  value = curr->value;
                }   

		if (curr->element_type == 'V') {
			vector_e[id_counter + 1 + voltageSources] += value;
			voltageSources++;
		}
		else if (curr->element_type == 'L') { 
		   	voltageSources++;
                }   
		else if (curr->element_type == 'I') {
			if (curr->terminalMNA_a != -1)			
				vector_e[curr->terminalMNA_a] -= value;
			if (curr->terminalMNA_b != -1)			
				vector_e[curr->terminalMNA_b] += value;
		}
                
        }

        matrixVectorMultiplication(C_with_h, vector_x, vector_b);

        for (i = 0; i < matrix_dim; i++) {
          vector_b[i] = vector_e[i] + vector_b[i];
        } 

}




//This function calculates the matrix A of the no sparse transient linear system, using Trapezodial method.

void TR_no_sparseA(double h) {

	int i, j;
	int matrix_dim = id_counter + 1 + g2_elements; 
	
	for (i = 0; i < matrix_dim; i++) {
		for (j = 0; j < matrix_dim; j++) {
			C_with_h[i][j] = C_no_sparse[i][j] * (2.0 / h);   
			matrix_A[i][j] = G_no_sparse[i][j] + C_with_h[i][j];
			C_with_h[i][j] = C_no_sparse[i][j] * (-2.0 / h);     
			C_with_h[i][j] = -(G_no_sparse[i][j] + C_with_h[i][j]);  
		}
	}     
	
	for (i = 0; i < matrix_dim; i++) {
		vector_e[i] = vector_b_dupl[i];
	}    
	
}



//This function calculates vector b of the no sparse transient linear system, using Trapezodial method

void TR_no_sparse_b(double time, double final_time) {

	int matrix_dim = id_counter + 1 + g2_elements;
        int voltageSources = 0;
 
	struct element *curr;

	int i;
	float value;   

        for (i = 0; i < matrix_dim; i++) {
          vector_e_prev[i] = vector_e[i]; 
          vector_e[i] = 0.0; 
        }

        // calculate e(t) which changes at every step of transient analysis  
	for (curr = root; curr != NULL; curr = curr->nxt_element) {

		if (curr->trans_src != NULL) {

			if (curr->trans_src->type == 1) {
				value = EXP_func(curr->trans_src->trans_structs_union.trans_EXP.i1, curr->trans_src->trans_structs_union.trans_EXP.i2,
                                                 curr->trans_src->trans_structs_union.trans_EXP.td1, curr->trans_src->trans_structs_union.trans_EXP.tc1,
                                                 curr->trans_src->trans_structs_union.trans_EXP.td2, curr->trans_src->trans_structs_union.trans_EXP.tc2, 
                                                 time, final_time); 
                        }              
 			else if (curr->trans_src->type == 2) { 
				value = SIN_func(curr->trans_src->trans_structs_union.trans_SIN.i1, curr->trans_src->trans_structs_union.trans_SIN.ia,
                                                 curr->trans_src->trans_structs_union.trans_SIN.fr, curr->trans_src->trans_structs_union.trans_SIN.td,
                                                 curr->trans_src->trans_structs_union.trans_SIN.df, curr->trans_src->trans_structs_union.trans_SIN.ph,
                                                 time, final_time);
                        }
                        else if (curr->trans_src->type == 3) { 
				value = PULSE_func(curr->trans_src->trans_structs_union.trans_PULSE.i1, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.i2, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.td, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.tr, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.tf, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.pw, 
                                                   curr->trans_src->trans_structs_union.trans_PULSE.per, 
                                                   time, final_time);

                        }
                        else if (curr->trans_src->type == 4) {
 				value = PWL_func(curr->trans_src->trans_structs_union.trans_PWL.n, 
                                                 curr->trans_src->trans_structs_union.trans_PWL.t, 
                                                 curr->trans_src->trans_structs_union.trans_PWL.i, 
                                                 time, final_time);
                        }
                        
		}   
		else {
                  value = curr->value;
                }   

		if (curr->element_type == 'V') {
			vector_e[id_counter + 1 + voltageSources] += value;
			voltageSources++;
		}
		else if (curr->element_type == 'L') { 
		   	voltageSources++;
                }   
		else if (curr->element_type == 'I') {
			if (curr->terminalMNA_a != -1)			
				vector_e[curr->terminalMNA_a] -= value;
			if (curr->terminalMNA_b != -1)			
				vector_e[curr->terminalMNA_b] += value;
		}
                
        }

        matrixVectorMultiplication(C_with_h, vector_x, vector_b); 

        for (i = 0; i < matrix_dim; i++) {
          vector_b[i] = vector_e_prev[i] + vector_b[i];
          vector_b[i] = vector_e[i] + vector_b[i];  
        }

}


void trans_free_no_sparse() {

        int i; 
	int matrix_dim = id_counter + 1 + g2_elements;

        for (i = 0; i < matrix_dim; i++) {
          free(G_no_sparse[i]); 
          free(C_no_sparse[i]);         
          free(C_with_h[i]); 
          free(matrix_A[i]);
        } 

        free(G_no_sparse);  
        free(C_no_sparse);
        free(C_with_h);
        free(matrix_A);  

	
	free(vector_e);
        free(vector_e_prev);
    	free(vector_b_dupl);


}
