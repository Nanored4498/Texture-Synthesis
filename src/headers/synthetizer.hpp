#ifndef SYNTHETIZER_DEF
#define SYNTHETIZER_DEF

#include "utils.hpp"

/**** Loading and initialisations ****/
/*
* The first line of parameters of each function contains parameters that will
* be used to compute parameters in the second line
* E is the brut sample image
* Ed is the downsized sample image if E was to big
* E2 is the torrified version of Ed if Ed was not a tore
* m, md and m2 are respective sizes of these three images
* new_E is true if E was created with new in C++ or with malloc in C
* ...
*/
int load_image(const char* filename, double to_tor,
				/** ----> **/
				uchar *&E, int &m, double &is_tore, double &new_E, uchar *&Ed, int &md);
void init_variables(uchar *Ed, int md, bool tor, const char* file,
					/** ----> **/
					int &m2, uchar*&E2, bool &have_folder, char* folder, int &L);
void init_live(int W0, int H0, uchar* E2, int m2, int L,
				/** ----> **/
				Pix *S[], int W[], int H[], uchar* El[]);
void init_live_WH(int L, int W0, int H0,
					/** ----> **/
					Pix *S[], int W[], int H[]);

// Algortihm for texture syntesis
// Return a map towars the pixel of E
Pix* synthesize (
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

// Apply the step l of the algorithm
void synthesize_step (
	int l,									// The step you want to apply
	Pix *S[], int W[], int H[],				// The array of maps and arrays of width and height of generated image at each step
	uchar *El[],							// The gaussian stack
	int md, int m2,							// The size of the (optionnaly downsized) sample and the size of the possibly torrified sample
	VD &r,									// The vector of amplitudes for the jitter step
	int L,									// Number of steps needed (log2 m)
	bool have_folder, char* folder,			// Some parameters usefull to find coherence files
	bool compute_co, int c, double kappa,	// Same parameters as in the synthetize function
	bool saveE = true						// If true then we write the image E in the file "out.png" else we write the map S
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