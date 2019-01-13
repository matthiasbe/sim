#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void read_matrix(char *filename, int size[2], double **A) {
	char buffer[1024] ;
	char *record,*line;
	int i=0,j=0;
	size[0] =0;
	size[1] =0;

	FILE *fstream = fopen(filename,"r");
	if(fstream == NULL)
	{
		printf("\n file opening failed ");
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

