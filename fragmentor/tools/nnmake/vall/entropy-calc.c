#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_LINE_LEN 225


/*
 * entropy-calc.c
 * Program for calculating profile entropy values for each residue in a Vall
 * database.
 * 
 * James Thompson <tex@u.washington.edu>
 * Copyright 2006 
 */

int aa_to_int( char aa );
char int_to_aa( int index );

// frequency matrix for amino acids. Order goes like this:
// A C D E F G H I K L M N P Q R S T V W Y
unsigned long aa_counts[20];
double ln_aa_freqs[20]; // stores the natural log of the aa frequencies
char* aa_order = "ACDEFGHIKLMNPQRSTVWY";
int num_aa = 20;
long total_aa = 0;

int main( int argc, char **argv ) {
    // process command-line arguments
    if ( argc != 2  ) {
        fprintf(stderr,"usage: entropy-calc vallfile\n");
        exit(1);
    }

    char* vallfile = argv[1];
    FILE *file = 0;
  
    // open the file read-only
    file = fopen( (const char*) vallfile,"r" );
    if ( file == NULL ) {
         fprintf(stderr,"Error opening file %s!\n", vallfile );
         exit(1);
    }

    char* line  = malloc( sizeof(char) * MAX_LINE_LEN  );
    
    // calculate amino acid frequencies
    if ( fgets(line,MAX_LINE_LEN,file) != NULL ) {
        char aa = line[6];
        aa_counts[aa_to_int(aa)]++;
    } else {
        fprintf(stderr,"Error reading from file %s!\n", vallfile );
        exit(1);
    }
    
    while ( fgets(line,MAX_LINE_LEN,file) != NULL ) {
        char aa = line[6];
        aa_counts[aa_to_int(aa)]++;
        total_aa++;
    }

    // fseek back to the beginning of the file, calculate entropies for each
    // profile and print them out.
    
    // E = sum( Pi * ln( Pi / Qi ) )
    // Where E is entropy, Qi is the background frequency of this
    // amino acid in the entire Vall, and Pi is the frequency of this
    // amino acid in this position.

    int i;
    double ln_total = log(total_aa);
    for ( i = 0; i < num_aa; i++ ) {
        ln_aa_freqs[i] = log(aa_counts[i]) - ln_total;
    }

    fseek(file, 0, SEEK_SET);
    while (  fgets(line,MAX_LINE_LEN,file) != NULL ) {
        double entropy = 0;

        // start at column 105
        int start = 104, size = 6;
        char* temp = malloc( sizeof(char) * size );
        for ( i = 0; i < num_aa; i++ ) {
            strncpy( temp, &line[start + i * size], size);
            double pi = atof(temp);
            entropy  += pi * ( pi - ln_aa_freqs[i] );
        }

        entropy = entropy / num_aa;

        printf("%f\n", entropy);
    }

    // free up temporary variables, close file handle
    free( line );
    fclose( file );

    return 0;
}

int aa_to_int( char aa ) {
    int i;
    for ( i = 0; i < num_aa; i++ ) {
        if ( aa == aa_order[i] ) {
            return i;
        }
    }

    return -1;
}

char int_to_aa( int index ) {
    return aa_order[index];
}
