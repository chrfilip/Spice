#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "parser.h"
#include "mna.h"



int main(int argc, char *argv[]) {

	double time;
	
	if(argc < 2){
		printf("Usage: ./parser filename\n");
		exit(1);
	}

	//Parse input file
	if ( parser(argv[1]) ) 
		printf("\nParsing OK!\n-------------------------------------\n");
	else {
		printf("Parsing error!\n");
		exit(1);
	}
	
	//Map node names to integers
	if(node_maping_hash())
		printf("\nNode maping OK!\n-------------------------------------\n");
	else {
		printf("Node maping error!\n");
		exit(1);
	}
	
	printf("number of nodes -> %d\nnumber of g2 elements -> %d\n", id_counter, g2_elements);
		
	printf("-------------------------------------\n");
	//Map plot sources
	if(plot)
		map_plot_nodes();
	
	if(sparse_option) {
		
		printf("Solve Ax=b system for DC analysis using sparse methods\n-------------------------------------\n");
		
		if(nz_approximation())
			printf("Nonzero approximation completed\n-------------------------------------\n");
		else {
			printf("Nonzero approximation failed\n");
			exit(1);
		}
		
		//Solve Ax=b system for DC analysis using sparse methods
		if(SPD_option) {
			if(ITER_option) {
				printf("Solving Ax=b using CG!\n-------------------------------------\n");
				if(DC_sweep){
				  DC_sweep_solver();
				}
				else{
				  if(sparse_CG(compression_function(create_triplet()))){
					printf("Sparse CG calculation completed\n");
				  }
				  else {
					printf("Sparse CG calculation failed\n"); 
					exit(1);
				  }
				}
				
				printf("CG algorithm completed\n-------------------------------------\n"); 
			}
			else {
				printf("Solving Ax=b using Cholesky!\n-------------------------------------\n");
				if(DC_sweep){
				  DC_sweep_solver();
				}
				else{
				  if(sparse_cholesky(compression_function(create_triplet()))){
					printf("Sparse Cholesky calculation completed\n");
				  }
				  else {
					printf("Sparse Cholesky calculation failed\n"); 
					exit(1);
				  }
				}
				
				
				printf("Cholesky algorithm completed\n-------------------------------------\n"); 
			}	
		}
		else if( !SPD_option && ITER_option) {
			printf("Solving Ax=b using Bi-CG!\n-------------------------------------\n");
			if(DC_sweep){
			  DC_sweep_solver();
			}
			else{
			  if(sparse_BiCG(compression_function(create_triplet()))){
				printf("Sparse Bi-CG calculation completed\n");
			  }
			  else {
				printf("Sparse Bi-CG calculation failed\n"); 
				exit(1);
			  }
			}
			
			
			printf("Bi-CG algorithm completed\n-------------------------------------\n"); 
		}
		else {
			printf("Solving Ax=b using LU!\n-------------------------------------\n");
			if(DC_sweep){
			  DC_sweep_solver();
			}
			else{
			  if(sparse_LU(compression_function(create_triplet()))){
				printf("Sparse LU calculation completed\n");
			  }
			  else {
				printf("Sparse LU calculation failed\n"); 
				exit(1);
			  }
			}
			
			printf("\tLU algorithm completed\n-------------------------------------\n");
		} 
		
	}
	else {
		
		printf("Solve Ax=b system for DC analysis using direct methods\n"); 
		
		//Circuit MNA analysis
		if( mna_analysis() )
			printf("MNA analysis OK!\n-------------------------------------\n");
		else {
			printf("MNA analysis Error! \n");
			exit (1);
		}
		//Solve Ax=b system for DC analysis
		if(SPD_option) {
			if(ITER_option) {
				printf("Solving Ax=b using CG!\n");
				if(DC_sweep){
				  DC_sweep_solver();
				}
				else{
				  if(cgSolver()){
					printf("CG algorithm completed!\n");
				    
				  }
				  else {
					printf("CG algorithm failed!\n");
				  }
				}
			}
			else {
			      printf("Solving Ax=b using Cholesky!\n");
			      cholesky_factorization();
			      if(DC_sweep){
				  DC_sweep_solver();
			      }
			      else{
				linearSolver();
			      }

			}	
		}
		else if( !SPD_option && ITER_option) {
			printf("Solving Ax=b using Bi-CG\n");
			if(DC_sweep){
			  DC_sweep_solver();
			}
			else{
			  if(BiCGSolver()){
				printf("Bi-CG algorithm completed\n");
			  }
			  else{
				printf("Bi-CG algorithm failed\n");
			  }
			}
			
		}
		else {
			printf("Solving Ax=b using LU!\n-------------------------------------\n");
			LU_factorization();
			if(DC_sweep) 
				DC_sweep_solver();
			else	
				linearSolver();
			
			printf("LU algorithm completed\n-------------------------------------\n");
		}
	}
	
	if (trans_option) {
		printf("sparse opt: %d\n", sparse_option);
			trans_init_no_sparse(); 
			
			GC_calculation();
			
			printf("C and G calculation for no sparse transient analysis completed\n");
			
			if (trans_method == BEuler) {
				
				printf("Calculation using 'Backward_Euler' for no sparse transient analysis started\n");
				
				BE_no_sparseA(tran_flags.time_step);
				
				if ( LU_factorization() ) {
					printf("LU factorization for no sparse transient analysis completed\n");
				}
				else {
					printf("LU factorization for no sparse transient analysis failed\n");
					exit(1); 
				}
				
				tran_output = fopen("tran_output.txt", "w");
				fclose(tran_output);
				
				for (time = 0.0; time < tran_flags.fin_time; time += tran_flags.time_step) {
					BE_no_sparse_b(time, tran_flags.fin_time);
					tran_output = fopen("tran_output.txt", "a");
					fprintf(tran_output, "Time: %lf\n", time);
					fclose(tran_output);
					linearSolver();

				}
				if(plot && trans_option){
				  plot = trans_option = 0;
				  linearSolver();
				  
				}
				
			}
			else if (trans_method == Trapezoidal) {
				
				printf("Calculation using 'Trapezoidal' for no sparse transient analysis started\n");
				
				TR_no_sparseA(tran_flags.time_step);
				
				if ( LU_factorization() ) {
					printf("LU factorization for no sparse transient analysis completed\n");
				}
				else {
					printf("LU factorization for no sparse transient analysis failed\n");
					exit(1); 
				}
				
				tran_output = fopen("tran_output.txt", "w");
				fclose(tran_output);
				
				for (time = 0.0; time < tran_flags.fin_time; time += tran_flags.time_step) {
					TR_no_sparse_b(time, tran_flags.fin_time);
					tran_output = fopen("tran_output.txt", "a");
					fprintf(tran_output, "Time: %lf\n", time);
					fclose(tran_output);
					linearSolver();
					
				}
				if(plot && trans_option){
				  plot = trans_option = 0;
				  linearSolver();
				  
				}
				
			}
			
			trans_free_no_sparse();
			
		
	}
	
	return 0;

}
