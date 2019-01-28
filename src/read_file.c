#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/** Read matrix from file in the simpler format :
 * the i-e line is the i-e line of the matrix, 
 * and each element of the line is space-separated
 */
void read_matrix(char *filename, int size[2], double **A) {
	char buffer[1024] ;
	char *record,*line;
	int i=0,j=0;
	size[0] =0;
	size[1] =0;

	// First, determine the matrix size
	FILE *fstream = fopen(filename,"r");
	if(fstream == NULL)
	{
		printf("file opening failed\n");
		exit(-1);
	}
	while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
	{
		j = 0;
		record = strtok(line," ");
		while(record != NULL)
		{
			j++;
			record = strtok(NULL," ");
		}
		++i ;
		if (j > size[1]) size[1] = j;
	}


	size[0] = i;
	i = 0;

	printf("detected matrix : %d x %d\n", size[0], size[1]);

	// Then fill mat with matrix element
	double (*mat)[size[1]] = (double (*)[]) malloc(sizeof(double)*size[0]*size[1]);
	fstream = fopen(filename,"r");
	while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
	{
		j = 0;
		record = strtok(line," ");
		while(record != NULL)
		{
			mat[i][j++] = atof(record);
			record = strtok(NULL," ");
		}
		++i ;
	}
	*A = mat[0];
}

/** Read a matrix from file, using the Matric Market format (.mtx)
 * representing sparse boolean matrices
 */
void read_mtx(char *filename, int size[2], double **A) {
	char buffer[1024];
	char *record, *line;

	FILE *fstream = fopen(filename, "r");
	if(fstream == NULL)
	{
		printf("file opening failed\n");
		exit(-1);
	}

	// Ignore first lines beginning with %
	do {
		line = fgets(buffer,sizeof(buffer),fstream);
	} while(buffer[0] == '%');

	// The next describe the size of the matrix
	size[0] = atoi(strtok(line, " "));
	size[1] = atoi(strtok(NULL, " "));

	printf("detected matrix : %d x %d\n", size[0], size[1]);

	int i=0,j=0;

	// Initialize all element to zero
	double (*mat)[size[1]] = (double (*)[]) malloc(sizeof(double)*size[0]*size[1]);
	for (int k = 0; k<size[0]; k++) {
		for (int l = 0; l<size[1]; l++) {
			mat[k][l] = 0;
		}
	}

	// Set to 1 referenced elements : a line gives an element's coords
	while((line=fgets(buffer,sizeof(buffer),fstream))!=NULL)
	{
		// First coord is the line
		i = atoi(strtok(line, " "));
		// Second coord is the column
		j = atoi(strtok(NULL, " "));
		mat[i][j] = 1;
	}
	*A = mat[0];
}

