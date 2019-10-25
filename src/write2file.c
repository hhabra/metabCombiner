#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

/*
 * File Description:
 * ------------------
 *
 *
 *
*/

/*
* Function: write2file
* -----------------------
*
* Takes a list of character-separated lines and prints to a file.
* 
* PARAMETERS:
*
* lines: character vector containing lines to be printed
*
* file: character string describing file destination
*
* groups: feature m/z group placement
* 
*/
void write2file(SEXP lines, SEXP file, SEXP groups)
{
	int *rgroups = INTEGER(groups);
	int n = LENGTH(groups);

	const char *filename = CHAR(STRING_ELT(file, 0));
	
	FILE *output = fopen(filename, "a");
	
	for(int i = 1; i < n; i++){
		if(rgroups[i] != rgroups[i-1])
			fputs("\n", output);
			
		const char *line = CHAR(STRING_ELT(lines, i));
	
		fputs(line, output);
		fputs("\n", output); 
	}

	fclose(output);
	
	return;
} 

