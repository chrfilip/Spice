#define STRING_SIZE 32
//#define PRINT_OUT

#define PWL_UPPER_LIMIT 10

//Transient methods
#define Trapezoidal 1
#define BEuler 0

struct tran_flags_t{
	double time_step;
	double fin_time;
};

struct tran_flags_t tran_flags;

// union of structs for all the types of transient sources
struct trans_src_t{

	int type;
	
	union {
		struct trans_EXP_t {
			double i1;
			double i2;
			double td1;
			double tc1;
			double td2;
			double tc2;
		} trans_EXP;
		
		struct trans_SIN_t {
			double i1;
			double ia;
			double fr;
			double td;
			double df;
			double ph;
		} trans_SIN;

		struct trans_PULSE_t {
			double i1;
			double i2;
			double td;
			double tr;
			double tf;
			double pw;
			double per;
		} trans_PULSE;

		struct trans_PWL_t{
			int n;
			double t[PWL_UPPER_LIMIT];
			double i[PWL_UPPER_LIMIT];
		} trans_PWL;

	} trans_structs_union;

};

//This data structure describes a seperate circuit element
struct element {

	char element_type;	//Element type
       	char *name;	//[STRING_SIZE];		//element name
	char *terminal_a;		//<node.+>
	char *terminal_b;		//<node.->
	char *terminal_c;		//MOS and BJT
	float value;		
	double length, width;	//Used in MOS
	float area;		//Used in BJT and Diode

	int terminalMNA_a;	//Used in MNA analysis
	int terminalMNA_b;

	struct trans_src_t *trans_src;

	struct element *nxt_element;
};


struct plot_node {
	
	char original_name[STRING_SIZE];
	int mapped_id;
	struct plot_node *next;
};

struct plot_node *plot_root;

//Options
bool SPD_option;
bool ITER_option;
bool DC_sweep;
bool sparse_option;
bool trans_option;			// turns true if a .TRAN statement exists
bool trans_method;			//0 for BE, 1 for TR
bool plot;					//true if plot statement



char dc_source[STRING_SIZE],source[STRING_SIZE];
float dc_start, dc_end, dc_step;
FILE *dc_output;
FILE *tran_output;
int pn;
char *temp_plot;
float itol;
int numberOfElements;

struct element *root;		//Head of the linked list

//Parsing functions
int parser(char *file);
void add_element(char elem_type, char *elem_name, char *fieldA, char *fieldB, char *fieldC, float elem_value, float elem_area, double elem_length, double elem_width, struct trans_src_t *transient_src);
void print_elements();
