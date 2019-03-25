#ifndef UTILS_DEF
#define UTILS_DEF

#include <vector>

/*******************/
/*** Definitions ***/
/*******************/

typedef unsigned char uchar;
typedef std::vector<double> VD;
typedef std::vector<int> VI;

struct Pix {
	int x, y;
	Pix(int a=0, int b=0): x(a), y(b) {}
	inline Pix operator+(const Pix &a) { return {x+a.x, y+a.y}; }
	inline Pix operator-(const Pix &a) { return {x-a.x, y-a.y}; }
	inline Pix operator%(int d) { return {x % d, y % d}; }
	inline Pix operator/(int d) { return {x / d, y / d}; }
};
inline Pix operator*(int a, const Pix &p) { return {a*p.x, a*p.y}; }
typedef std::vector<Pix> VP;
typedef std::vector<VP> VVP;

/********************/
/*** Neighborhood ***/
/********************/

// fill N with the k-neighborhood of (x, y) in the map (S -> E) of size W*H
// where E is a squared image of size m.
void getNeighb_S(int x, int y, int k, Pix* S, int W, int H, uchar* E, int m, VI &N);
// fill N with the k-neighborhood of (x, y) in the squared image E of size m
// with footsteps h
void getNeighb(int x, int y, int k, uchar* E, int m, VI &N, int h);
// return the distance between two k-neighborhood a and b
inline int distNeighb(VI &a, VI &b, int k) {
	int k2 = k*k*3;
	int res = 0, diff;
	for(int i = 0; i < k2; i++) {
		diff = a[i] - b[i];
		res += diff * diff;
	}
	return res;
}

/***********************/
/*** Shape of images ***/
/***********************/

// return true iff the squared image im of size m is a tore
bool is_tor(uchar* im, int m);
// return a squared image of size 2m which the image im in a corner
// with miror images of im in other corners to obtain a tore
uchar* torrify(uchar* im, int m);
// return a squared image of size min(W, H) from an image im of size W*H
uchar* square(uchar* im, int W, int H);

/*********************/
/*** Saving images ***/
/*********************/

// Save the map (S -> E) of size W*H where E is a squared image of size m
// in the file name
void save(Pix* S, int W, int H, uchar* E, int m, const char* name);
// Same as the function above but apply a filter on borders of S
void save_smooth(Pix* S, int W, int H, uchar* E, int m, const char* name);
// Save the map S of size W*H as a 2D image red/green in the file name
void saveS(Pix* S, int W, int H, int m, const char* name);

#endif