#ifndef WAVReader_h
#define WAVReader_h

#include <iostream>
#include <fstream>
#include <cmath>

#define BUFF_SIZE 1024
#define META_SIZE 44
#define NUM_HARMS 5
#define _USE_MATH_DEFINES

class WAVReader
{
  public:
    //file - file to read
    //trim - number of bytes to skip at the beginning of file
    //numIts - number of averages to make per sample period
    //fund - fundamental frequecy
    //numHarms - number of harmonics to store
    //numSamples - max amount of data points
    //output - array of arrays of data points
    //	each data point is an array of doubles, standing for a0,b0,
    //	a1,b1,a2,b2....a(numHarms - 1),b(numHarms - 1), where a and b
    //	are the fourier coefficients
    //returns the number of data points read
    int read(const char * file, int trim, int numIts, double f0, 
    	int maxSamples, double * outp);
    //normalizes data so that the total amplitude is 1
    void normalize(double * data, int dataPts);
};


#endif
