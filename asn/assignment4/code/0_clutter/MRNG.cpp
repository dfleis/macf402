#include "MRNG.h"
#include <iostream>
#include <math.h> // pow

#define NORM 2.328306549295728e-10
#define M1   4294967087.0
#define M2   4294944443.0
#define A12     1403580.0
#define A13N     810728.0
#define A21      527612.0
#define A23N    1370589.0

double MRNG::runif(double s[]) {
    long k;
    double p1, p2;
	
    p1 = A12 * s[1] - A13N * s[0];
	k = p1 / M1;
	p1 -= k * M1;
	if (p1 < 0.0) p1 += M1;
	s[0] = s[1]; s[1] = s[2]; s[2] = p1;

	p2 = A21 * s[5] - A23N * s[3];
	k  = p2 / M2;
	p2 -= k * M2;
	if (p2 < 0.0) p2 += M2;
	s[3] = s[4]; s[4] = s[5]; s[5] = p2;

    double r;
    
	if (p1 <= p2) r = ((p1 - p2 + M1) * NORM);
	else r = ((p1 - p2) * NORM);
    return r;
}













