#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>

#define MAX_SIZE 30000
#define PI 3.1415962
#define FACTOR 1
#define ITS 1

using namespace std;



// convert two bytes to one double in the range -1 to 1
double bytesToDouble(char firstByte, char secondByte)
{
    // convert two bytes to one short (little endian)
    short s = (short)((secondByte << 8) | firstByte);
    // convert to range from -1 to (just below) 1
    return s / 32768.0;//s > 0 ? 1 - s / 32768.0 : -s / 32768.0 - 1;
}


// Returns left and right double arrays. 'right' will be null if sound is mono.
void openWav(double * left, double * right, int maxSize)
{
		ifstream myFile;// (filename, ios::in | ios::binary);
		myFile.open("violin_A4.wav");
		char wav[100000];
		myFile.read(wav, 100000);
		int read = myFile.gcount();
		
    // Determine if mono or stereo
    int channels = wav[22];     // Forget byte 23 as 99.999% of WAVs are 1 or 2 channels


    // Get past all the other sub chunks to get to the data subchunk:
    int pos = 12;   // First Subchunk ID from 12 to 16


    // Keep iterating until we find the data chunk (i.e. 64 61 74 61 ...... (i.e. 100 97 116 97 in decimal))
    while (!(wav[pos] == 100 && wav[pos + 1] == 97 && wav[pos + 2] == 116 && wav[pos + 3] == 97))
    {
        pos += 4;
        int chunkSize = wav[pos] + wav[pos + 1] * 256 + wav[pos + 2] * 65536 + wav[pos + 3] * 16777216;
        pos += 4 + chunkSize;
    }
    pos += 8;


    // Pos is now positioned to start of actual sound data.
    int samples = (read - pos) / 2;     // 2 bytes per sample (16 bit sound mono)
    if (channels == 2) samples /= 2;        // 4 bytes per sample (16 bit stereo)


    // Allocate memory (right will be null if only mono sound)
    //left = new double[samples];
    //if (channels == 2) right = new double[samples];
    //else right = null;

		//skip first 3 kb bytes
		
		pos += 1024 * 10;
		
    // Write to double array/s:
    int i = 0;
    while (i < maxSize)
    {
        * left = bytesToDouble(wav[pos], wav[pos + 1]);
        left++;
        pos += 2;
        if (channels == 2)
        {
            * right = bytesToDouble(wav[pos], wav[pos + 1]);
            right++;
            pos += 2;
        }
        i++;
    }
}


void printWave(double * wv, int num, int inc)
{
	for (int i = 0; i < num; i++)
	{
		int amp = 20;
		for (int l = 0; l < amp + min((* wv) * amp * FACTOR, 0.0); l++)
		{
			cout << " ";
		}
		for (int l = amp + min((* wv) * amp * FACTOR, 0.0); l < 
			amp + max((* wv) * amp * FACTOR, 0.0); l++)
		{
			cout << "-";
		}
		cout << endl;
		wv += inc;
	}
}

double an(double * wv, double T0, double dt, int n)
{
	double sm = 0;
	for (int i = 0; i <= T0 / dt * ITS; i++)
	{
		sm += 2 / T0 * dt * (* wv) * FACTOR * cos(n * i * dt / T0 * 2 * PI);
		wv++;
	}
	return sm / ITS;
}

double bn(double * wv, double T0, double dt, int n)
{
	double sm = 0;
	//cout << "--------------" << endl;
	for (int i = 0; i <= T0 / dt * ITS; i++)
	{
		//cout << "wv : " << (* wv) << " i : " << i << " x : " <<  i * dt / T0 * 2 * PI << endl;
		//cout << "sin sq: " <<  << endl;
		sm += 2 / T0 * dt * (* wv) * FACTOR * sin(n * i * dt / T0 * 2 * PI);
		//cout << "wv: " << (* wv) <<
		//	" sin(nx) " << sin(n * i * dt / T0 * 2 * PI) <<
		//	" wv*sin(nx)*dt " << dt * (* wv) * sin(n * i * dt / T0 * 2 * PI) <<
		//	" sum: " << sm << endl;
		wv++;
	}
	return sm / ITS;
}



int main()
{
	double leftAr[MAX_SIZE];
	double rightAr[MAX_SIZE];
	openWav(leftAr, rightAr, MAX_SIZE);
	
	double * left = leftAr;
	double * right = rightAr;
	
	//cout << leftAr[0] << " " << leftAr[1] << " " << leftAr[2] << endl;
	
	//printWave(left, 45, 1);
	
	
	//testing an algorithm
	//double T0 = 2 * PI;//1.0 / 1000;//fundamental period
	//double dt = 2 * PI / 20;//1.0 / 2000 / 24;//time step of data point in file
	
	//double newAr[24];
	//for (int i = 0; i < 24; i++)
	//{
	//	newAr[i] = cos(i * 4 * PI / 20);
	//}
	
	double T0 = 1.0 / 440;//fundamental period
	double dt = 1.0 / 22050;//time step of data point in file
	
	double thresh = .03;
	
	string graph[20];
	double step = .001;
	for (int i = 0; i < 20; i++)
	{
		graph[i] = "";
	}
	
	
	for (int i = 0; i < 20; i++) 
	{
		
		double aCoef = an(left, T0, dt, i);
		//if (aCoef * aCoef > thresh * thresh) 
		//{
		//cout << "a" << i << " = " << aCoef << endl;
		//}
		double bCoef = bn(left, T0, dt, i);
		//if (bCoef * bCoef > thresh * thresh)
		//{
		//cout << "b" << i << " = " << bCoef << endl;
		//}
		
		double cnt = 0;
		for (int j = 0; j < 20; j++)
		{
			graph[19 - j] += aCoef * aCoef >= cnt * cnt ? "|" : " "; 
			graph[19 - j] += bCoef * bCoef >= cnt * cnt ? "|" : " ";
			cnt += step;
		}
	}
	for (int i = 0; i < 20; i++)
	{
		cout << graph[i] << endl;
	}
	//cout << left[0] << endl << left[1] << endl << left[2] << endl;
	return 0;
}
