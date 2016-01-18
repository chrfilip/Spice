#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "parser.h"


int parser(char *file) {
  
	FILE *f;
	char line[256];
	char temp_plot[32];
	char option[8];
	char name[STRING_SIZE], terminal_a[STRING_SIZE], terminal_b[STRING_SIZE], terminal_c[STRING_SIZE];
	char elem_type;
	int i=0, line_counter = 0;
	float value, area = 1;
	double width, length;
	char tokens[5][10];
	char delims[] = " =\n";
	char *result = NULL;
	int k=0;
	int temp_length = 0;
	struct trans_src_t *transient_src;

	pn = 0;
	//Initialize trans pointer
	transient_src = NULL;

	SPD_option = 0;
	DC_sweep = 0;
	ITER_option = 0;
	sparse_option = 0;
	trans_option = 0;
	trans_method = Trapezoidal; //default transient method
	plot = 0;
	
	plot_root = NULL;
	root = NULL;
	
	itol = 0.001; // default value
	numberOfElements = 0;

	//Open netlist file for reading
	f = fopen(file,"r");
	if(f == NULL){
		printf("ERROR: Could not open netlist!\n");
		return 0;
	}

	struct plot_node *new_plot_node = NULL;
	
	//Get next line
	while( fgets(line, sizeof(line), f) != NULL)	//read a line
	{
		line_counter++;		//keeps line number
		
		if(strlen(line) <= 2)		//Is the line empty?
		{
			continue;
		}
		
		
		//Remove accidental spaces at the beggining of the line	
		if(line[i]== ' ')
		{
			while(line[i] == ' ') 
				i++;
		}

		if(line[i] == '\n')		//Check if there is a comment
		{
			continue;
		}
		if(line[i] == '\r')		//Check if there is a comment
		{
			continue;
		}
		if(line[i] == '*')		//Check if there is a comment
		{
			continue;
		}

		//Lowercase to uppercase
		if(line[i] >= 'a' && line[i] <= 'z')
			line[i] = line[i] - 'a' +'A';

		//Check for . statements
		if(line[i] == '.') {
			
		  
			//Check for Transient statements
			if(line[i+1] == 'T') {
				
				if( sscanf( &line[i+5], " %lf %lf", &(tran_flags.time_step), &(tran_flags.fin_time)) != 2 )
				{
					printf("ERROR: Incomplete description in line %d!\n", line_counter);
					return 0;
				}
				//printf("trans option %lf, %lf\n", tran_flags.time_step, tran_flags.fin_time);
				
				trans_option = 1;
				continue;
			}
			else if(line[i+1] == 'O') {
				
				
				if(strncmp(&line[i+9], "METHOD", 6)==0) {
				
					if(strncmp(&line[i+16], "BE", 2)==0)
						trans_method = BEuler;
					else if(strncmp(&line[i+16], "TR", 2)==0)
						trans_method = Trapezoidal;
					//printf("Method ok! = %d\n", trans_method);
					continue;
				}
				
				result = strtok( line, delims );
				while( result != NULL ) {
					strcpy(tokens[k], result);
					//printf( "result is %s\n", tokens[k] );
					result = strtok( NULL, delims );
					k++;
				}
				
				switch(k){
					
					case 2:
						if(strcmp(tokens[1],"SPD")==0)
							SPD_option = 1;
						else if(strncmp(tokens[1],"SPARSE",6)==0)
							sparse_option = 1;
						else if(strcmp(tokens[1],"ITER")==0)
							ITER_option = 1;
						break;
						
					case 3:
						if(strcmp(tokens[1],"SPD")==0 && strcmp(tokens[2],"ITER")==0) {
							SPD_option = 1;
							ITER_option = 1;
						}
						else if(strcmp(tokens[1],"SPARSE")==0 && strcmp(tokens[2],"SPD")==0) {
							sparse_option = 1;
							SPD_option = 1;
						}
						else if(strcmp(tokens[1],"SPARSE")==0 && strcmp(tokens[2],"ITER")==0) {
							sparse_option = 1;
							ITER_option = 1;
						}
						break;
					case 4: 
						if(strcmp(tokens[1],"ITER")==0 && strcmp(tokens[2],"ITOL")==0) {
							ITER_option = 1;
							itol = strtod(tokens[3], NULL);
							printf("Itol: %f\n", itol);  
						}
						
						else if(strcmp(tokens[1],"SPARSE")==0 && strcmp(tokens[2],"SPD")==0 && strcmp(tokens[3],"ITER")==0) {
							sparse_option = 1;
							ITER_option = 1;
							SPD_option = 1;
						}
						break;
					case 5:	
						if (strcmp(tokens[1],"SPD")==0 && strcmp(tokens[2],"ITER")==0 && strcmp(tokens[3],"ITOL")==0) {
							ITER_option = 1;
							SPD_option = 1;
							itol = strtod(tokens[4],NULL);
							printf("Itol: %f\n", itol);
						}
						else if(strcmp(tokens[1],"SPARSE")==0 && strcmp(tokens[2],"ITER")==0 && strcmp(tokens[3],"ITOL")==0) {
							sparse_option = 1;
							ITER_option = 1;
							itol = strtod(tokens[4],NULL);
							printf("Itol: %f\n", itol);
						}
						break;
					case 6: 
						if(strcmp(tokens[1],"SPARSE")==0 && strcmp(tokens[2],"SPD")==0 && strcmp(tokens[3],"ITER")==0 && strcmp(tokens[4],"ITOL")==0) {
							sparse_option = 1;
							SPD_option = 1;
							ITER_option =1;
							itol = strtod(tokens[5],NULL);
                                                        printf("Itol: %f\n", itol);   
						}
						break;
				}

			} 
			else if( line[i+1] == 'D') {
				
				if(sscanf(line, "%s %s %f %f %f", option, dc_source, &dc_start, &dc_end, &dc_step)==1){
				  DC_sweep=0;
				}
				else if(sscanf(line, "%s %s %f %f %f", option, dc_source, &dc_start, &dc_end, &dc_step)!=5) {
					printf("ERROR: Incomplete description in line %d!\n", line_counter);
					return 0;
				}			
				if(strcmp(option,".DC") == 0 && sscanf(line, "%s %s %f %f %f", option, dc_source, &dc_start, &dc_end, &dc_step)==5){
					DC_sweep = 1;
				}
			}
			else if (line[i+1] == 'P'){
				if (sscanf(line,"%s %s", option, source)!= 2){
					printf("ERROR: Incomplete description in line %d!\n", line_counter);
					return 0;
				}
				if (DC_sweep == 0 && trans_option == 0){
					printf("WARNING: Missing DC sweep or TRAN option, plotting failed!\n");
					continue;
				}
				
				
				while(line[i++]!=' ');
				
				if(trans_option || DC_sweep) {
					
					while( i < strlen(line) && sscanf(&line[i], "%*[V]%*[(]%[^)]%*[)]%*[ ,\t,\n]", temp_plot)) {
					
						//printf("temp plot: %s\n", temp_plot);
						new_plot_node = (struct plot_node *)malloc(sizeof(struct plot_node));
						strcpy(new_plot_node->original_name, temp_plot);
						//printf("plot: %s\n", new_plot_node->original_name);
						
						new_plot_node->next = plot_root;
						plot_root = new_plot_node;
						pn++;
						i+= strlen(new_plot_node->original_name)+4;
					}
					
					plot = 1;
					
				}

			}
			continue;
		}
		else if( line[i] == 'V' || line[i] == 'I' || line[i] == 'R' || line[i] == 'C' || line[i] == 'L')
		{	
			elem_type = line[i];
			line[i] = ' ';
			if( sscanf( line, "%s %s %s %e", name, terminal_a, terminal_b, &value) != 4 )
			{
				printf("ERROR: Incomplete description in line %d!\n", line_counter);
				return 0;
			}
	
			if( strlen(name) == 1 && name[0] == '0' ){
			
				printf("ERROR: Incorrect element name for line: %c%s %s %s %.2e \n", elem_type, name, terminal_a, terminal_b, value);
				printf("Element name '0' is reserved for ground\n");
				return 0;
			}

			if( terminal_a == terminal_b){

				printf("The element %c%s is not connected properly (short-circuit)\n", elem_type, name);
				return 0;
			}
 
			if( elem_type == 'V' || elem_type == 'I') {
				
				temp_length = (int)strlen(name) + (int)strlen(terminal_a) + (int)strlen(terminal_b) + 4;

				while(line[temp_length++]!=' ');

				if(temp_length < strlen(line) && strncmp(&line[temp_length], "EXP", 3) == 0){

									
					transient_src = (struct trans_src_t*)malloc(sizeof(struct trans_src_t));
					if(transient_src == NULL) {
						printf("Memory allocation error in line %d\n", line_counter);
						return 0;
					}
					
					transient_src -> type = 1;
					if( sscanf(&line[temp_length+3], " (%lf %lf %lf %lf %lf %lf)", &transient_src->trans_structs_union.trans_EXP.i1, &transient_src->trans_structs_union.trans_EXP.i2, &transient_src->trans_structs_union.trans_EXP.td1, &transient_src->trans_structs_union.trans_EXP.tc1, &transient_src->trans_structs_union.trans_EXP.td2, &transient_src->trans_structs_union.trans_EXP.tc2) != 6) {
						
						printf("ERROR: Incomplete description in line %d! (EXP) \n", line_counter);
						return 0;
					}

					
				}
				
				else if(temp_length < strlen(line) && strncmp(&line[temp_length], "SIN", 3) == 0){

					transient_src = (struct trans_src_t*)malloc(sizeof(struct trans_src_t));
					if(transient_src == NULL) {
						printf("Memory allocation error in line %d\n", line_counter);
						return 0;
					}
					transient_src -> type = 2;
					if( sscanf(&line[temp_length+3], " (%lf %lf %lf %lf %lf %lf)", &transient_src->trans_structs_union.trans_SIN.i1, &transient_src->trans_structs_union.trans_SIN.ia, &transient_src->trans_structs_union.trans_SIN.fr, &transient_src->trans_structs_union.trans_SIN.td, &transient_src->trans_structs_union.trans_SIN.df, &transient_src->trans_structs_union.trans_SIN.ph) != 6) {
						
						printf("ERROR: Incomplete description in line %d! (SIN) \n", line_counter);
						return 0;
					}


				}
				
				else if(temp_length < strlen(line) && strncmp(&line[temp_length], "PULSE", 5) == 0){

					transient_src = (struct trans_src_t*)malloc(sizeof(struct trans_src_t));
					if(transient_src == NULL) {
						printf("Memory allocation error in line %d\n", line_counter);
						return 0;
					}
					transient_src -> type = 3;
					if( sscanf(&line[temp_length+5], " (%lf %lf %lf %lf %lf %lf %lf)", &transient_src->trans_structs_union.trans_PULSE.i1, &transient_src->trans_structs_union.trans_PULSE.i2, &transient_src->trans_structs_union.trans_PULSE.td, &transient_src->trans_structs_union.trans_PULSE.tr, &transient_src->trans_structs_union.trans_PULSE.tf, &transient_src->trans_structs_union.trans_PULSE.pw, &transient_src->trans_structs_union.trans_PULSE.per) != 7) {
						
						printf("ERROR: Incomplete description in line %d! (PULSE) \n", line_counter);
						return 0;
					}

					
				}
				
				else if(temp_length < strlen(line) && strncmp(&line[temp_length], "PWL", 3) == 0){

					transient_src = (struct trans_src_t*)malloc(sizeof(struct trans_src_t));
					if(transient_src == NULL) {
						printf("Memory allocation error in line %d\n", line_counter);
						return 0;
					}
					transient_src -> type = 4;
					transient_src->trans_structs_union.trans_PWL.n = 0;
					temp_length += 3;
					
					while(temp_length < strlen(line)) {

						if(sscanf(&line[temp_length], " (%lf %lf)", &transient_src->trans_structs_union.trans_PWL.t[transient_src->trans_structs_union.trans_PWL.n], &transient_src->trans_structs_union.trans_PWL.i[transient_src->trans_structs_union.trans_PWL.n]) != 2) {
						
							printf("ERROR: Incomplete description in line %d!\n", line_counter);
							return 0;
						}
						transient_src->trans_structs_union.trans_PWL.n++;
						temp_length++;
						while(line[temp_length++] != ' ');
						while(line[temp_length++] != ' ');
						temp_length--;
						
						if(transient_src->trans_structs_union.trans_PWL.n == PWL_UPPER_LIMIT) {
							printf("Maximum number of PWL values! \n");
							return 0;
						}
					}

				}
				else {
					transient_src = NULL;
				}
				
				
			}
			strcpy(terminal_c, "NULL");
			length = width = -1;
			area = -1;
		}
		else if( line[i] == 'Q' )
		{
			elem_type = line[i];
			line[i] = ' ';
			if( sscanf( line, "%s %s %s %s %f", name, terminal_a, terminal_b, terminal_c, &area ) != 5 )
			{
				printf("ERROR: Incomplete description in line %d!\n", line_counter);
				return 0;
			}

			if( terminal_a == terminal_b && terminal_b == terminal_c ){

				printf("The element %c%s is not connected properly (short-circuit)\n", elem_type, name);
				return 0;
			}
			value = -1;
			length = width = -1;
		
		}
		else if( line[i] == 'M' )
		{
			elem_type = line[i];
			line[i] = ' ';
			if( sscanf( line, "%s %s %s %s %lf %lf", name, terminal_a, terminal_b, terminal_c, &length, &width ) != 6 )
			{
				printf("ERROR: Incomplete description in line %d!\n", line_counter);
				return 0;
			}

			if( terminal_a == terminal_b && terminal_b == terminal_c ){

				printf("The element %c%s is not connected properly (short-circuit)\n", elem_type, name);
				return 0;
			}
			value = area = -1;

		}
		else if(line[i] == 'D')
		{
			elem_type = line[i];
			line[i] = ' ';
			if( sscanf( line, "%s %s %s %f", name, terminal_a, terminal_b, &area ) != 4 )
			{
				printf("ERROR: Incomplete description in line %d!\n", line_counter);
				return 0;
			}

			if( terminal_a == terminal_b ){

				printf("The element %c%s is not connected properly (short-circuit)\n", elem_type, name);
				return 0;
			}
			strcpy(terminal_c, "NULL");
			length = width = -1;
			value = -1;

		}
		else
		{
			printf("ERROR: Unknown element, line %d\n", line_counter);
			return 0;
		}

		//Add element in the list
		add_element( elem_type, name, terminal_a, terminal_b, terminal_c, value, area, length, width, transient_src);
		
		transient_src = NULL;
		
		numberOfElements++;
		i=0;

	}

	return 1;

}


void add_element(char elem_type, char *elem_name, char *fieldA, char *fieldB, char *fieldC, float elem_value, float elem_area, double elem_length, double elem_width, struct trans_src_t *trans_src)
{


	struct element *current;

	current = (struct element *)malloc(sizeof(struct element));

	current->element_type = elem_type;
	current->name = strdup(elem_name);
	current->terminal_a = strdup(fieldA);
	current->terminal_b = strdup(fieldB);
	current->terminal_c = strdup(fieldC);
	current->value = elem_value;
	current->area = elem_area;
	current->length = elem_length;
	current->width = elem_width;
	current->terminalMNA_a = -1;
	current->terminalMNA_b = -1;
	current->trans_src = trans_src;

	
	current->nxt_element = root;
	root = current;

}


void print_elements(){

	struct element *current;

	for(current = root; current != NULL; current= current->nxt_element){

		if(current->element_type == 'V' || current->element_type == 'R' || current->element_type == 'I' || current->element_type == 'C'  || current->element_type == 'L'){
			
			printf("%c%s %s %s %.2e\n", current->element_type, current->name, current->terminal_a, current->terminal_b, current->value);

		}
		else if( current->element_type == 'Q'){

			printf("%c%s %s %s %s %f \n", current->element_type, current->name, current->terminal_a, current->terminal_b, current->terminal_c, current->area);
			
		}
		else if( current->element_type == 'M'){

			printf("%c%s %s %s %s %lf %lf\n", current->element_type, current->name, current->terminal_a, current->terminal_b, current->terminal_c, current->length, current->width );
		
		}
		else if( current->element_type == 'D'){

			printf("%c%s %s %s %f\n", current->element_type, current->name, current->terminal_a, current->terminal_b, current->area);

		}
	}

	printf("Number of elements: %d\n", numberOfElements);

}
