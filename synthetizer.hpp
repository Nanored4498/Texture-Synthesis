#ifndef SYNTHETIZER_DEF
#define SYNTHETIZER_DEF

#include "utils.hpp"

// Algortihm for texture syntesis
// Return a map to the pixel of E
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

// Reduce the size of im by dividing the width of the image by 2^p
// m is the original size of im which is a squared image
uchar* downsample(uchar* im, int m, int p);
// Increase the quality of an image
uchar* magnify(
	int ml,		// The size of the downsampled image used to create the map S
	uchar* Eh,	// The original sample image before the downsizing
	int mh,		// The original width of the sample image
	Pix* S,		// The map to the downsized sample image
	int W,		// The width of the map S
	int H,		// The height of the map S
	int &Wh,	// The new width will be stored in this integer
	int &Hh		// The new height will be stored in this integer
);

#endif