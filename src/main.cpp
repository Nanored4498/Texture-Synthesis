#include <iostream>
#include <cstring>
#include <cmath>
#include "stb_image.h"
#include "stb_image_write.h"
#include "utils.hpp"
#include "synthetizer.hpp"

int main(int argc, char* argv[]) {
	/****************************/
	/******** Parameters ********/
	/****************************/
	int W = 2, H = 2;
	int c = 2;
	double kappa = 0.08;
	// VD r = {0.2, 0.3, 0.32, 0.3, 0.2, 0.1, 0.1, 0.1, 0.0};
	VD r = {0.65, 0.35, 0.48, 0.35, 0.3, 0.2, 0.2, 0.1, 0.0};
	bool compute_co = false;
	bool to_tor = false;
	/****************************/
	/****************************/

	// Parsing
	if(argc < 2) {
		std::cerr << "This is a texture synthetiser." << std::endl;
		std::cerr << "You can use it by typing:" << std::endl;
		std::cerr << "\t" << argv[0] << " [filename] [-c] [-t]" << std::endl;
		std::cerr << "Where:" << std::endl;
		std::cerr << " -filename is the name of the image file used as an example" << std::endl;
		std::cerr << " -c if this parameter is given then coherence is computed if it doesn't already exist" << std::endl;
		std::cerr << " -t if this parameter is given then the image in input is torified" << std::endl;
		return 1;
	}
	for(int i = 2; i < argc; i++) {
		if(strcmp(argv[i], "-c") == 0) compute_co = true;
		else if(strcmp(argv[i], "-t") == 0) to_tor = true;
	}
	char* filename = argv[1];

	uchar *E=0, *Ed=0;
	int m, md;
	double is_tore, new_E;
	if(load_image(filename, to_tor, E, m, is_tore, new_E, Ed, md))
		return 1;

	// Check if we have to downsize example
	if(md < m) {
		Pix* S = synthesize(Ed, md, r, c, kappa, W, H, !is_tore, filename, compute_co);
		int Wh, Hh;
		uchar* Sh = magnify(md, E, m, S, W, H, Wh, Hh);
		stbi_write_png("magnific.png", Wh, Hh, 3, Sh, 0);
		delete[] S;
		delete[] Sh;
		delete[] Ed;
	} else {
		Pix* S = synthesize(E, m, r, c, kappa, W, H, !is_tore, filename, compute_co);
		int L = 1 << (int) ceil(log2(m));
		save_smooth(S, W*L, H*L, E, m, "magnific.png");
		delete[] S;
	}

	if(new_E) delete[] E;
	else free(E);
	return 0;
}