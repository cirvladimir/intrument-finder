#include "WAVReader.h"

#define NUM_DATS 1
#define NUM_ROWS 25

using namespace std;

void graph(double * dt, string tit)
{
	for (int i = 0; i < NUM_HARMS * 3 - tit.size() / 2 - 1; i++)
	{
		cout << "-";
	}
	cout << tit;
	for (int i = 0; i < NUM_HARMS * 3 - (tit.size() + 1) / 2; i++)
	{
		cout << "-";
	}
	cout << endl;
	double step = 1.0 / NUM_ROWS;
	double cur = 1;
	for (int i = 0; i < NUM_ROWS; i++)
	{
		for (int j = 0; j < NUM_HARMS * 2; j++)
		{
			if (abs(dt[j]) >= cur)
			{
				if (abs(dt[j]) - step < cur)
					cout << "__ ";
				else
					cout << "|| "; 
			}
			else
				cout << "   ";
		}
		cout << endl;
		cur -= step;
	}
	for (int i = 0; i < NUM_HARMS; i++)
	{
		cout << "a" << i << " b" << i << " ";
	}
	cout << endl;
}

void print(double * dat)
{
	for (int i = 0; i < NUM_HARMS; i++)
	{
	 	cout << "a" << i << " = " << dat[i * 2] << endl;
	 	cout << "b" << i << " = " << dat[i * 2 + 1] << endl; 
	}
}

int main()
{
	cout << "Welcome to sound grapher" << endl;
	WAVReader r;
	double dat[NUM_DATS * NUM_HARMS * 2];
	const char * insts[4] = { "singer_A4.wav", "flute_A4.wav", 
		"piano_A4.wav", "violin_A4.wav" };
	for (int i = 2; i < 3; i++)
	{
		int numPts = r.read(insts[i], 5000, 1, 550, NUM_DATS, dat);
		r.normalize(dat, numPts);
		graph(dat, insts[i]);
		//print(dat);
	}
}
