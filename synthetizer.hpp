#ifndef SYNTHETIZER_DEF
#define SYNTHETIZER_DEF

#include "utils.hpp"

Pix* synthesize(
	uchar* E,					// The sample image which have to be a square
	int m,						// The size of the image
	VD &r,						// The vector of amplitudes for the jitter step
	int c,						// The number of correction passes
	double kappa,				// The penalty for the coherence usage
	int W,						// The width of the result (this is the number of duplications of E so the real width is W*m)
	int H,						// The height of the result (this is the number of duplications of E so the real height is H*m)
	bool tor=false,				// If the sample E has been torified
	const char* file="",		// The folder in which are or will be stored coherence images
	bool compute_co=false		// If true and coherence images don't exist then compute them
);

uchar* downsample(uchar* im, int m, int p);
uchar* magnify(int ml, uchar* Eh, int mh, Pix* S, int W, int H);

#endif