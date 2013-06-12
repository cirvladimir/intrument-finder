#include "WAVReader.h"

using namespace std;

int readLEnd(char * beg, int num)
{
	int outp = 0;
	for (int i = 0; i < num; i++)
	{
		outp += (*(beg++)) << (i * 8);
	}
	return outp;
}

//file - file to read
//trim - number of bytes to skip at the beginning of file
//numIts - number of averages to make per sample period
//fund - fundamental frequecy
//numSamples - max amount of data points
//output - array of arrays of data points
//	each data point is an array of doubles, standing for a0,b0,
//	a1,b1,a2,b2....a(numHarms - 1),b(numHarms - 1), where a and b
//	are the fourier coefficients
//returns the number of data points read
int WAVReader::read(const char * file, int trim, int numIts, double f0, 
	int maxSamples, double * outp)
{
	ifstream myFile;
	myFile.open(file);
	
	//meta data:
	char meta[META_SIZE];
	myFile.read(meta, META_SIZE);
	
	int numChanels = meta[22] | (meta[23] << 8);
	int sampleRate = readLEnd(&meta[24], 4);
	int sampBitSize = readLEnd(&meta[34], 2);
	
	int sampByteSize = readLEnd(meta + 32, 2);
	
	//trimming beginning bytes
	myFile.seekg(trim * sampByteSize, ios_base::cur);
	
	//reading samples
	int datsRead = 0;
	int bRead;
	char wav[BUFF_SIZE];
	
	//number of samples to read to use for each period * numIts
	double numPerDat = 1 / f0 * sampleRate * numIts;
	double sampsLeft = numPerDat;
	int sampsRead = 0;
	
	double harms[NUM_HARMS * 2];
	for (int i = 0; i < NUM_HARMS * 2; i++)
	{
		harms[i] = 0;
	}
		
	do
	{
		myFile.read(wav, BUFF_SIZE);
		bRead = myFile.gcount();
		
		char * endPos = wav + bRead;
		char * curPos = wav;
		
		//loop through all bytes read in buffer
		while ((curPos < endPos) && (datsRead < maxSamples))
		{
			//read value of sound, increment pointer
			short locVal;
			if (sampBitSize == 16)
			{
				locVal = * curPos;
				locVal = locVal | (*(curPos + 1) << 8);
			}
			else
			{
				locVal = (* curPos << 8);
				locVal /= 256;
			}
			curPos += sampByteSize;
			
			//loop through all harmonics, add to fourier coefficient
			for (int n = 0; n < NUM_HARMS; n++)
			{
				//adding to an
				harms[2 * n] += 2 * f0 / numIts * locVal * 
					cos(n * 2 * M_PI * f0 * sampsRead / sampleRate) / sampleRate;
			
				//adding to bn
				harms[2 * n + 1] += 2 * f0 / numIts * locVal * 
					sin(n * 2 * M_PI * f0 * sampsRead / sampleRate) / sampleRate;
			}
		
			sampsRead++;
			sampsLeft--;
			//if finished sampling period * numIts, write that to
			//	output, reset counter
			if (sampsLeft < 0) 
			{
				//writing to outp
				for (int i = 0; i < NUM_HARMS * 2; i++)
				{
					*(outp++) = harms[i];
				}
				for (int i = 0; i < NUM_HARMS * 2; i++)
				{
					harms[i] = 0;
				}
				//resetting counter
				sampsLeft += numPerDat;
				datsRead++;
			}
		}
	} while ((datsRead < maxSamples) && (bRead == BUFF_SIZE));
	return datsRead;
}

void WAVReader::normalize(double * data, int dataPts)
{
	for (int i = 0; i < dataPts; i++)
	{
		double sum = 0;
		for (int j = 0; j < NUM_HARMS * 2; j++)
		{
			sum += abs(data[j]);
		}
		for (int j = 0; j < NUM_HARMS * 2; j++)
		{
			*(data) = (*data) / sum;
			data++;
		}
	}
}
