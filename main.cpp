#include <iostream>
#include <cstring>
#include <cmath>
#include "stb_image.h"
#include "stb_image_write.h"
#include "utils.hpp"
#include "synthetizer.hpp"

int main(int argc, char* argv[]) {

	// Parameters
	int W = 3, H = 2;
	VD r = {0.2, 0.3, 0.35, 0.3, 0.2, 0.1, 0.1, 0.1};
	bool compute_co = false;
	bool to_tor = false;

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

	// Loading image
	char* filename = argv[1];
	int w, h, col;
	uchar* E = stbi_load(filename, &w, &h, &col, 3);
	if (!E) {
		std::cerr << "loading failed" << std::endl;
		return 1;
	}

	// Eventually resize to a square
	bool new_E = false;
	int m = w;
	if(w != h) {
		uchar* temp = square(E, w, h);
		m = std::min(w, h);
		free(E);
		E = temp;
		new_E = true;
	}

	// If the user want to torify image then torify the image
	if(to_tor) {
		uchar* tE = torrify(E, m);
		if(new_E) delete[] E;
		else free(E);
		E = tE;
		m *= 2;
	}

	// Check if the example is a tore
	bool is_tore = is_tor(E, m);
	if(is_tore)
		std::cerr << "The example has been considered as a tore" << std::endl;

	// Check if we have to downsize example
	if(m > 256) {
		int ml = m;
		int ds = 0;
		while(ml > 256) ml >>= 1, ds++;
		uchar* d = downsample(E, m, ds);
		Pix* S = synthesize(d, ml, r, 3, 0.2, W, H, !is_tore, filename, compute_co);
		int Wh, Hh;
		uchar* Sh = magnify(ml, E, m, S, W, H, Wh, Hh);
		stbi_write_png("magnific.png", Wh, Hh, 3, Sh, 0);
		delete[] S;
		delete[] Sh;
		delete[] d;
	} else {
		Pix* S = synthesize(E, m, r, 3, 0.2, W, H, !is_tore, filename, compute_co);
		int L = 1 << (int) ceil(log2(m));
		save_smooth(S, W*L, H*L, E, m, "magnific.png");
		delete[] S;
	}

	if(new_E) delete[] E;
	else free(E);
	return 0;
}