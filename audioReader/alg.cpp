#include "WAVReader.h"
#include "../linearMapping/data.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

//number of vectors to make per instrument
#define NUM_DATS 30
//number of sample instruments
#define NUM_INSTS 4
//number of periods to go through
#define NUM_PERIODS 15
//fundamental frequency
#define FUNDAMENTAL 550

int main()
{
	cout << "Welcome to sound analyzer" << endl;
	WAVReader r;
	double dat[NUM_DATS * NUM_HARMS * 2];
	
	//sample instruments
	const char * insts[NUM_INSTS] = { "singer_A4.wav", "flute_A4.wav", 
		"piano_A4.wav", "violin_A4.wav" };
	
	//all data points stored here
	double * allPts = new double[NUM_DATS * NUM_INSTS * NUM_HARMS * 2];
	//used to store next data to allPts without changing the pointer
	double * allPtsWriter = allPts;
	
	//average instruments
	double ** avgPts = new double*[NUM_INSTS];
	//setting average insruments to zero vectors
	for (int i = 0; i < NUM_INSTS; i++)
	{
		avgPts[i] = new double[NUM_HARMS * 2];
		for (int j = 0; j < NUM_HARMS * 2; j++)
		{
			avgPts[i][j] = 0;
		}
	}
	
	//gathering data
	int totPts = 0;
	for (int instInd = 0; instInd < NUM_INSTS; instInd++)
	{
		cout << "Instrument " << instInd << " - " << insts[instInd] << endl;
		//read some points
		//store to central point list
		//record average for instrument
		int numPts = r.read(insts[instInd], 5000, NUM_PERIODS, FUNDAMENTAL, NUM_DATS, dat);
		totPts += numPts;
		r.normalize(dat, numPts);
		
		for (int dtPt = 0; dtPt < numPts * NUM_HARMS * 2; dtPt++)
		{
			avgPts[instInd][dtPt % (NUM_HARMS * 2)] += dat[dtPt] / numPts;
			*(allPtsWriter++) = dat[dtPt];
		}
	}
	
	//load data/compute eigen values
	Data data(NUM_HARMS * 2);
	data.load(allPts, totPts);
	
	//compute orthornormal basis
	data.compOrthoBasis(.1);
	cout << endl;
	
	//control experiment: try matching sample files used with instruments
	cout << "Control: testing against same files" << endl;
	int testContStarts[3] = { 5000, 8000, 15000 };
	for (int contTestInd = 0;  contTestInd < 3; contTestInd++) 
	{
		cout << "testing starting at sample " << testContStarts[contTestInd] << endl;
		for (int i = 0; i < 4; i++)
		{
			//read/normalize a vector of fourier coefficients
			double sampData[NUM_HARMS * 2];
			r.normalize(sampData, r.read(insts[i], testContStarts[contTestInd], NUM_PERIODS, FUNDAMENTAL, 1, sampData));
			
			//match vector to one of the average instrument vectors
			int match = data.match(avgPts, 4, sampData);
			cout << "instument " << i << " matched with instrument " 
				<< match << endl;
		}
	}
	cout << endl;
	
	//experiment: try matching another violin file with instruments used in sample
	const char * tFiles[3] = {"Violin-vib-in-A4.wav", "mystery2.wav", "mystery1.wav"};
	for (int i = 0; i < 3; i++)
	{
		const char * tFile = tFiles[i];
		cout << "Testing " << tFile  << endl;
	
		int testActStarts[4] = { 5000, 15000, 25000, 55000 };
		for (int testInd = 0; testInd < 4; testInd++)
		{
			cout << "Testing starting at sample " << testActStarts[testInd] << endl;
			//read/normalize a vector of fourier coefficients
			double difSamp[NUM_HARMS * 2];
			r.normalize(difSamp, r.read(tFile, testActStarts[testInd], NUM_PERIODS, 550, 1, difSamp));
		
			//match vector with one of the average instrument vectors
			int match = data.match(avgPts, 4, difSamp);
			cout << "Unknown instument matched with instrument " 
				<< match << endl;
		}
	}
}
