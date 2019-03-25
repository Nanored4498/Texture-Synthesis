#ifndef COHERENCE_DEF
#define COHERENCE_DEF

#include "utils.hpp"

/*****************************************/
/**** Computation of coherence images ****/
/*****************************************/

// For each pixel of E in the square of size m, add a new nearest neighbor in C.
// E is an image of size m2 which can be different from m because E can be a
// torified image of the original sample image.
void appendCoherence(VVP &C, uchar* E, int m, int m2, int h);

// Initialise the coherence C with identity in the square of size m
void initCoherence(VVP &C, int m);

// Load the coherence strored in the file filename
// and add the nearest neighbors in C
void loadCoherence(VVP &C, const char* filename);

// Write the kth nearest neighbors stored in the matrix C in the image
// filename. Coherence is stored in a square of size m
void writeCoherence(VVP &C, int k, const char* filename, int m);

#endif