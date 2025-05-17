// 210101062_digitrecognition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <sstream>
#include <Windows.h>
using namespace std;
#pragma comment(lib, "winmm.lib")

//constants
const int framesize = 320;
const int limit = 5000;
const int P = 12;
const double PI = 3.142857142857;
const int CB_SIZE = 32;				// Desired codebook size
const double EPSILON = 0.03;		// Perturbation factor
const int M = 32;		// No.of obsevation symbols per state
const int N = 5;		// No.of states
const int MAX_T = 160;
const int TRAIN_FILES = 30;

//global variables
double dcShift = 0;
double nFactor = 1;
int num_samples = 0;
int file_no = 0;
double silenceEnergy = 0;
int multiplier = 100;
long double P_O_given_lambda = 0;
int steadyStart=0;
int steadyEnd=0;

double tokhuraWeights[]={1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
double samples[100000] = {0};
double energy[10000] = {0};
double steady_frame[100000]={0};
double R[P+1] = {0};
double Ai[P+1] = {0};
double C[P+1] = {0};
int file_frames[350]={0};
int O[MAX_T+1]={0};
int T=0;
int q_star[MAX_T+1];	//state sequence.
long double p_star = 0, prev_p_star = -1;
long double alpha[MAX_T+1][N+1];
long double beta[MAX_T+1][N+1];
long double gamma[MAX_T+1][N+1];
long double delta[MAX_T+1][N+1];
int Psi[MAX_T+1][N+1]; 
long double xi[MAX_T+1][N+1][N+1];

//Model parameters A, B and Pi
long double A[N+1][N+1], B[N+1][M+1], Pi[N+1];
long double Abar[N+1][N+1], Bbar[N+1][M+1], Pibar[N+1];
long double a_average[N+1][N+1], b_average[N+1][M+1], pi_average[N+1];
long double threshold = 1e-30;

FILE *uni=NULL; 
short int waveIn[16025 * 3];
vector<vector< double>> universe;
vector<vector<double>> codebook;

//function to reset the current model to the unbiased model
void resetModel(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j] = 0;
		}
	}
	A[1][1] = 0.8;
	A[1][2] = 0.2;
	A[2][2] = 0.8;
	A[2][3] = 0.2;
	A[3][3] = 0.8;
	A[3][4] = 0.2;
	A[4][4] = 0.8;
	A[4][5] = 0.2;
	A[5][5] = 1;
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			B[i][j] = 1.0/32*1.0;
		}
	}

	for(int i=1; i<=N; i++){
		Pi[i] = 0;
	}
	Pi[1]=1;
}

//function to reset the average model to zeroes
void resetAvgModel(){
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			a_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			b_average[i][j] = 0;
		}
	}

	for(int i=1; i<=N; i++){
		pi_average[i] = 0;
	}
}

//function to compute the DC offset of the speech signal
void getDCShift(){
	dcShift = 0;
	for(int i=0; i<num_samples; i++)
	{
		dcShift += samples[i];
	}
	dcShift /= num_samples;
}

//function to compute the normalization factor of the speech signal
void getNFactor(){
	nFactor = 1;
	double max = INT_MIN;
	for(int i=0; i<num_samples; i++)
	{
		if(abs(samples[i]) > max){
			max = abs(samples[i]);
		}
	}
	nFactor = (double)limit/max;
}

//function to compute the energy of a frame
double computeEnergy(double* arr,int frame_size){
    double energy = 0;
    for(int i=0; i<frame_size; i++) 
	{
        energy += (arr[i]*arr[i]);
        if (energy < 0) {
           printf("Overflow detected!");
           return -1;
        }
    }
	energy /= frame_size;
    return energy;
}

//function to calculate the silence threshold energy
void computeEnergyThresh(){
	FILE *f1=NULL;				//Read the data from a text file
	int err = fopen_s( &f1, "silence.txt", "r");
	if(err != NULL)
	{
		printf("\nCannot open the file\n");
		system("pause");
		exit(1);
	}                      
	
	double x=0;
	int silenceCnt=0;
	double silenceSamples[100000]={0};		//read the file
	while( !feof(f1) )
	{
		if( fscanf_s( f1, "%lf", &x) == 1){
			silenceSamples[ silenceCnt ] = x;
			silenceCnt++;
		}
		else{
			char line[1000];
			fgets(line,100,f1);
		}
	}

	int silenceFrames=0;
	silenceEnergy=0;
	for(int i=0; i<silenceCnt; i+=framesize)	//finding the avg energy
	{
		int frame_size = (i + framesize <= num_samples) ? framesize : num_samples - i;
		silenceEnergy += computeEnergy(silenceSamples+i,frame_size);
		silenceFrames++;
	}
	silenceEnergy/=silenceFrames;
}

//function to compute the Ri values
void computeRis(double* frame, int frame_size){
	for(int m=0; m<=P; m++){
		R[m]=0;
		for(int k=0; k<frame_size-m; k++){
			R[m] += frame[k]*frame[k+m];
		}
	}
}

//function to compute the Ai values using Levinson Durbin algorithm
void computeAis(){                    
	double Alpha[13][13],E[13],K[13];
	E[0] = R[0];				  //step 1

	for(int i=1; i<=P; i++){					
		double sum=0;
		for(int j=1; j<=i-1; j++){				
			sum += Alpha[i-1][j]*R[i-j];	
		}
			
		K[i]=(R[i]-sum)/E[i-1];		//step 2		
		Alpha[i][i]=K[i];			//step 3
		for(int j=1; j<=i-1; j++){
			Alpha[i][j]=Alpha[i-1][j] - K[i]*Alpha[i-1][i-j];		//step 4
		}
		E[i]=(1-(K[i]*K[i]))*E[i-1];		//step 5
	}
	
	for(int i=1;i<=P;i++){
		Ai[i]= Alpha[P][i];
	}
}

//function to calculate the cepstral coeff Ci's
void computeCis(){
	double sum=0;
	C[0]=log(R[0]*R[0]);

	for(int m=1; m<=P; m++){
		sum=0;
		for(int k=1; k<m; k++){
			sum += (k*C[k]*Ai[m-k])/(m*1.0);
		}
		C[m]=Ai[m]+sum;
	}
	//applying raised sin window on cis 
	/*for(int m=1;m<=P;m++){ 
		C[m]*=(P/2)*sin((PI*m)/P);
	}*/	
}

//funtion to find the steady frames
void findSteadyFrames(){
	int energySize = 0;
	double en = 0.0;

	for(int i=0; i<num_samples; i++) {
		//if the n is a multiple of N then we have to store energy to the array
		if(i%framesize == 0 && i!=0){
			en /= framesize;				//taking average
			energy[energySize++] = en;	//store the energy value
			en = 0;						//resetting en to 0
		}
		en += samples[i] * samples[i];
	}
	energy[energySize++] = en;

	int minSamples=9600;
	steadyStart = 0;
	steadyEnd = num_samples-1;
	for(int i=0; i<energySize-4; i++){
		if(steadyStart == -1 && steadyEnd == -1 && energy[i+1] > multiplier * silenceEnergy && energy[i+2] > multiplier * silenceEnergy && energy[i+3] > multiplier * silenceEnergy && energy[i+4] > multiplier * silenceEnergy){
			steadyStart = i*framesize;
		}else if(steadyStart != -1 && steadyEnd == -1 && energy[i+1] <= multiplier * silenceEnergy && energy[i+2] <= multiplier * silenceEnergy && energy[i+3] <= multiplier * silenceEnergy && energy[i+4] <= multiplier * silenceEnergy){
			int diff = i*framesize - steadyStart;
			if(diff < minSamples){
				steadyStart = 0 > (steadyStart - (minSamples - diff)/2) ? 0 : (steadyStart - (minSamples - diff)/2);
				steadyEnd = energySize*framesize < (i*framesize + (minSamples - diff)/2) ? energySize*framesize : (i*framesize + (minSamples - diff)/2);
			}
			else 
				steadyEnd = i*framesize;
		}else if(steadyStart != -1 && steadyEnd!= -1) break;
	}
}

//function to calculate the Tokhura distance between the current centroid and current row of universe
double tokhuraDistance(const vector<double> &currentCentroid, const vector<double> &currentVector) {
    double ans = 0.0;
    for (int i=0; i<P; i++) {
        double d = currentCentroid[i] - currentVector[i];
        ans += (tokhuraWeights[i] * d * d);			//using tokhura distance
    }
    return ans;
}

//function to calculate the centroid of the universe
vector<double> calculateCentroid(const vector<vector<double>>& universe) {
    vector<double> centroid(P, 0.0);
    for (int i=0;i<universe.size();i++) {
        for (int j=0; j<P; j++) {
            centroid[j] += universe[i][j];		//add all the values of universe
        }
    }
    for (int i=0; i<P; i++) {
        centroid[i] /= universe.size();			//take average for finding centroid
    }
    return centroid;
}

//function to double the size of the codebook by splitting each current code vector
void splitCodebook(vector<vector<double>>& codebook) {
    vector<vector<double>> newCodebook;
	vector <double> t1,t2;
	//cout<<codebook.size()<<endl;
    for (int j=0;j<codebook.size();j++) {
        //Apply perturbation for each vector component
        for (int i = 0; i <P; i++) {
			long double v1 = ((1 + EPSILON) * codebook[j][i]);
            t1.push_back(v1);
			long double v2 = ((1 - EPSILON) * codebook[j][i]);
            t2.push_back(v2);
        }

        //Add the perturbed vectors to the new codebook
        newCodebook.push_back(t1);
        newCodebook.push_back(t1);
    }
    codebook = newCodebook;  // Update the codebook with the newly split entries
}

//Lloyds algorithm to perform vector quantization
void kMeans(vector<vector<double>>& universe, vector<vector<double>>& codebook, double delta) {
    double prevDistortion = DBL_MAX;
    int m = 0;
    while(true) {
        vector<vector<vector<double>>> regions(CB_SIZE); // Regions for vectors
        double totalDistortion = 0.0;

        for (int j=0; j<universe.size(); j++) {			//Assign vectors to nearest codebook
            double minDist = DBL_MAX;
            int closestRegion = 0;

            for (int i = 0; i < codebook.size(); i++) {
                double dist = tokhuraDistance(universe[j], codebook[i]);
                if (dist < minDist) {
                    minDist = dist;
                    closestRegion = i;
                }
            }
            regions[closestRegion].push_back(universe[j]);
            totalDistortion += minDist;
        }
		double avgDistortion = totalDistortion / universe.size(); 
        
        for (int i = 0; i < codebook.size(); i++) {		//Calculate the new centroids for each region
            if (!regions[i].empty()) {
                codebook[i] = calculateCentroid(regions[i]);
            }
        }
		
        if (abs(prevDistortion - avgDistortion) < delta) {		//Check for convergence
            break;
        }
        prevDistortion = avgDistortion;
        m++;
		cout<<"Distortion in iteration "<<m<<": "<<prevDistortion<<endl;
    }
}

//function to generate observation sequence from the ci's
int genObservations(double* frame){
	int  min_index = 0;
	double min = DBL_MAX;
	double sum = 0 ;

	for (int j=0; j<CB_SIZE; j++){
		sum=0;
		for (int i=1; i<=P; i++){
			sum += tokhuraWeights[i] * (frame[i] - codebook[j][i-1])*(frame[i] - codebook[j][i-1]);
		}
		if (sum < min){
			min = sum;
			min_index = j;
		}
	}
	return (min_index + 1);
}

//function to calculate alpha variable to find the solution of problem number 1
void forwardProcedure(){
	long double sum;
	for(int i=1; i<=N; i++){			//initialization
		alpha[1][i] = Pi[i]*B[i][O[1]];
	}
	
	for (int t=1; t<T; t++){			//induction
		for (int j=1; j<=N; j++){
			sum = 0;
			for (int i=1; i<=N; i++){
				sum += alpha[t][i] * A[i][j];
			}
			alpha[t + 1][j] = sum * B[j][O[t + 1]];
		}
	}

	P_O_given_lambda = 0;
	for(int i=1; i<=N; i++){			//termination
		P_O_given_lambda += alpha[T][i];
	}	
}

//function to calculate beta variable
void backwardProcedure(){
	long double sum;
	for(int i=1; i<=N; i++){		//initialization
		beta[T][i] = 1.0;
	}
	for(int t=T-1; t>=1; t--){		//induction
		int obs = O[t+1];
		for(int i=1;i<=N;i++){
			sum = 0;
			for(int j=1; j<=N; j++){
				sum += B[j][obs]*A[i][j]*beta[t+1][j];
			}
			beta[t][i]=sum;
		}
	}
}

//function to find the single best state sequence for the given observation sequence
void viterbiAlgorithm(){
	p_star = 0;
    for(int i=1; i<=N; i++){			//Initialization
        delta[1][i] = Pi[i] * B[i][O[1]];
        Psi[1][i] = 0;
    }

	for(int j=1; j<=N; j++){			//Induction
		for(int t=2; t<=T; t++){
            long double maxVal = 0, temp = 0;
            int index = 0;            
            for(int i=1; i<=N; i++){
                temp = (delta[t-1][i] * A[i][j]);
                if(temp > maxVal){
					maxVal = temp;
					index = i;
				}
            }
            delta[t][j] = maxVal * B[j][O[t]];
			Psi[t][j] = index;
        }
    }

    long double max = 0;
    for(int i=1; i<=N; i++){			//Termination
        if(delta[T][i] > max) {
			max = delta[T][i];
			q_star[T] = i;
		}
        p_star = max;
    }

    for(int t=T-1; t>0; t--){			////Path (state sequence) backtracking
        q_star[t] = Psi[t+1][q_star[t+1]];
    }
}

//function to calculate xi variable
void computeXi(){
	long double denum;
	for(int t=1; t<=T-1; t++){
		denum = 0;
		for(int i=1; i<=N; i++){
			for(int j=1; j<=N; j++){
				denum += alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
			}
		}
		for(int i=1; i<=N; i++){
			long double z = 0;
			for(int j=1; j<=N; j++){
				z = alpha[t][i]*A[i][j]*B[j][O[t+1]]*beta[t+1][j];
				xi[t][i][j]= z/denum;
			}
		}
	}
}

//function to compute gamma values
void computeGamma(){
	for (int t=1; t<=T-1; t++){
		for (int i=1; i<=N; i++){
			gamma[t][i] = 0;
			for (int j=1; j<=N; j++){
				gamma[t][i] += xi[t][i][j];
			}
		}
	}
}

//function to find solution to the reestimation problem
void reestimateParameters(){
	long double sum1=0 , sum2 =0;
	for(int i=1; i<=N; i++){
		Pibar[i] = gamma[1][i];					//Reevaluating Pi
	}
	for(int i=1; i<=N; i++){					//Reevaluating A
		for(int j=1; j<=N; j++){
			sum1 = 0, sum2 = 0;
			for(int t=1; t<=T-1; t++){
				sum1 += xi[t][i][j];
				sum2 += gamma[t][i];
			}
			Abar[i][j] = sum1/sum2;
		}
	}
	for(int j=1; j<=N; j++){					//Reevaluating B
		for(int k=1; k<=M; k++){
			sum1 =0.0, sum2 =0.0;
			for(int t=1; t<T; t++){
				sum1 += gamma[t][j];
				if(O[t]==k){
					sum2 += gamma[t][j];				
				}
			}
			Bbar[j][k] = sum2/sum1;
		}

		//Applying threshold to B to prevent #IND problem
		int max_idx = 0;
		long double max_B = Bbar[j][1];
		for(int k=2; k<=M; k++){
			if(Bbar[j][k] > max_B){
				max_B = Bbar[j][k];
				max_idx = k;
			}
		}
		long double b_threshold = 0.0;
		for(int k=1; k<=M; k++){
			if(Bbar[j][k] <= threshold){
				b_threshold++;
				Bbar[j][k] += threshold;
			}
		}

		Bbar[j][max_idx] -= b_threshold *threshold;

	}

	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			A[i][j]= Abar[i][j];
			//cout<<A[i][j]<<" "; 
		}
		//cout<<endl;
		for(int k=1; k<=M; k++){
			B[i][k] = Bbar[i][k];
			//cout<<B[i][k]<<" ";
		}
		//cout<<endl;
		Pi[i] = Pibar[i];
		//cout<<Pi[i]<<endl;
	}
}

//function to add the current model values to avg model
void addModels(){
	int i, j;
	for (i=1; i<=N; i++){
		for (j=1; j<=N; j++){
			a_average[i][j] += A[i][j];
		}
	}
	for (i=1; i<=N; i++){
		pi_average[i] += Pi[i];
	}
	for (int i=1; i<=N; i++){
		for (int j=1; j<=M; j++){
			b_average[i][j] += B[i][j];
		}
	}
}

//function to compute average model
void computeAvgModel(int total_iterations){
	int i, j;
	for (i=1; i<=N; i++){
		for (j=1; j<=N; j++){
			a_average[i][j] /= total_iterations;

		}
	}
	for (i=1; i<=N; i++){
		for (j=1; j<=M; j++){
			b_average[i][j] /= total_iterations;
		}
	}
	for (i=1; i<=N; i++){
		pi_average[i] /= total_iterations;
	}
}

//function to write average models to text file
void writeAvgModels(int digit){
	char a_file_avg[100], b_file_avg[100], pi_file_avg[100], ind[3];

	sprintf(a_file_avg, "avgModels/digit_%d_A.txt", digit);
	FILE *fp = fopen(a_file_avg, "w");
	for(int i=1; i<=N; i++){
		for(int j=1; j<=N; j++){
			fprintf(fp, "%Le   ", a_average[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	
	sprintf(b_file_avg, "avgModels/digit_%d_B.txt", digit);
	ofstream fout(b_file_avg);
	for(int i=1; i<=N; i++){
		for(int j=1; j<=M; j++){
			fout<<b_average[i][j]<<"   ";
		}
		fout<<endl;
	}
	fout.close();

	sprintf(pi_file_avg, "avgModels/digit_%d_PI.txt", digit);
	fp = fopen(pi_file_avg, "w");
	for(int i=1; i<=N; i++){
		fprintf(fp, "%Le   ", pi_average[i]);
	}
	fclose(fp);
}

//function to read average models from the text file
void readAvgModels(int digit){
	char path[100];
	sprintf(path, "avgModels/digit_%d_A.txt", digit);
	fstream fin;
	fin.open(path);    
    
	if(!fin){			//file does not exist
		cout<<"Couldn't open file: "<<path<<"\n";
		return;
	}
	long double word;                     
	int row = 1, col = 1;
	while(fin >> word){			//read the file until input is available
        col = 1;
        A[row][col++] = word;
		//cout<<word<<" ~ ";
		//cout<<A[row][1]<<endl;
        for(int i=2; i<=N; i++){
            fin>>word;
            A[row][col++] = word;
        }
        row++;
    }
	fin.close();

	sprintf(path, "avgModels/digit_%d_B.txt", digit);
	fin.open(path);    

	if(!fin){				//file does not exist
		cout<<"Couldn't open file: "<<path<<"\n";
		return;
	}                    
	row = 1;
	col = 1;
	while(fin >> word){		//read the file
        col = 1;
        B[row][col++] = word;
        for(int i=2; i<=M; i++){
            fin>>word;
            B[row][col++] = word;
        }
        row++;
    }
	fin.close();

	sprintf(path, "avgModels/digit_%d_PI.txt", digit);
	fin.open(path);    

	if(!fin){				//file does not exist
		cout<<"Couldn't open file: "<<path<<"\n";
		return;
	}                    
	row = 1;
	while(fin >> word){		//read the file
        Pi[row++]=word;
        row++;
    }
	fin.close();
}

void writedataTofile(LPSTR lpData,DWORD dwBufferLength);

//function for testing a live recording
void testingLive(){
	cout<<"\n------------------------------------TESTING LIVE----------------------------------------\n";
	num_samples = 16025 * 3;
	for(int i=1; i<num_samples; i++){
		samples[i] = waveIn[i];
	}
	getDCShift();
	getNFactor();
	for(int i=0; i<num_samples; i++){
		samples[i]=((samples[i]-dcShift)*nFactor);
	}
	findSteadyFrames();
	int sampSize=0, k=steadyStart, T=0;
	while(k<steadyEnd){
		int frame_size = (k + framesize <= steadyEnd) ? framesize : steadyEnd - k;
		for(int h=0; h<frame_size; h++){
			steady_frame[sampSize++] = samples[k++];
		}
		computeRis(steady_frame+sampSize-frame_size, frame_size);
		computeAis();
		computeCis();
		O[++T] = genObservations(C);
	}
	file_frames[file_no]=T;
	int tesPred = 0;
	int flag = 0;
	double maxProbDigit = 0;
	//double pogivenLambda[10]={0};
	for(k = 0; k<=9; k++){
		readAvgModels(k);
		forwardProcedure();
		cout<<"P(O | model - "<<k<<") : "<<P_O_given_lambda<<endl;
		if(P_O_given_lambda > maxProbDigit){
			maxProbDigit= P_O_given_lambda;
			tesPred=k;
		}
		//pogivenLambda[k] = P_O_given_lambda;
		resetAvgModel();
	}
	cout<<"Predicted digit: "<<tesPred<<endl;
}

//function for testing the test recordings
void testRecordings(){
	char filename[100], line[100], test_file[100];
	int correctAns = 0, totalCorrectAns = 0;
	int file_no = 0;

	cout<<"\n------------------------------------TEST RECORDINGS----------------------------------------\n";
	for(int i=0;i<10;i++){
		correctAns=0;
		for(int j=31;j<=40;j++){
			char path[100];
			sprintf(path, "210101062_dataset/English/txt/210101062_E_%d_%d.txt", i, j);
			double x;
			int err;
			FILE *f1=NULL; // Read the data 
			cout<<path<<endl;
			err = fopen_s( &f1, (const char *)path, "r");
			if(err != NULL)
			{
				printf("\nCannot open the file\n");
				system("pause");
				exit(1);
			}                      
			num_samples=0;
			while( !feof(f1) )			//read the file
			{
				if( fscanf_s( f1, "%lf", &x) == 1)
				{
					samples[ num_samples ] = x;
					num_samples++;
				}
				else{
					//printf( "Error in line %d %d\n ", num_samples,file_no );
					char line[1000];
					fgets(line,100,f1);
				}
			}
			fclose(f1);
			getDCShift();
			getNFactor();
			for(int i=0;i<num_samples;i++){
				samples[i]=((samples[i]-dcShift)*nFactor);
			}
			findSteadyFrames();
			int sampSize=0,k=steadyStart,T=0;
			while(k<steadyEnd){
				int frame_size = (k + framesize <= steadyEnd) ? framesize : steadyEnd - k;
				for(int h=0; h<frame_size; h++){
					steady_frame[sampSize++] = samples[k++];
				}
				computeRis(steady_frame+sampSize-frame_size,frame_size);
				computeAis();
				computeCis();
				O[++T] = genObservations(C);
				//steady_frame[sampSize]*= 0.54-0.46*cos(2*PI*steady_frame[sampSize]/(framesize-1));
			}
			file_frames[file_no]=T;
			int testPred = 0;
			double maxProbDigit = 0;
			for(k = 0; k<=9; k++){
				readAvgModels(k);
				forwardProcedure();
				cout<<"P(O | model - "<<k<<") : "<<P_O_given_lambda<<endl;
				if(P_O_given_lambda > maxProbDigit){
					maxProbDigit= P_O_given_lambda;
					testPred = k;
				}
				resetAvgModel();
			}
			if(testPred == i) correctAns++, totalCorrectAns++;
			cout<<"Predicted digit: "<<testPred<<endl;
			file_no++;
			//cout<<"file "<<i<<" "<<j<<" completed\n";
		}
		printf("Accuracy for the digit %d is : %lf % \n", i, (correctAns / 10.0f)*100);
	}
	printf("Accuracy for the model is : %lf % \n",  (totalCorrectAns / 100.0f)*100);
}

void PlayRecord()
{
	 const int NUMPTS = 16025 * 3;   // 3 seconds
	 int sampleRate = 16025;  
	 // 'short int' is a 16-bit type; I request 16-bit samples below
	 // for 8-bit capture, you'd    use 'unsigned char' or 'BYTE' 8-bit types
	 HWAVEIN  hWaveIn;
	 WAVEFORMATEX pFormat;
	 pFormat.wFormatTag=WAVE_FORMAT_PCM;     // simple, uncompressed format
	 pFormat.nChannels=1;                    //  1=mono, 2=stereo
	 pFormat.nSamplesPerSec=sampleRate;      // 44100
	 pFormat.nAvgBytesPerSec=sampleRate*2;   // = nSamplesPerSec * n.Channels * wBitsPerSample/8
	 pFormat.nBlockAlign=2;                  // = n.Channels * wBitsPerSample/8
	 pFormat.wBitsPerSample=16;              //  16 for high quality, 8 for telephone-grade
	 pFormat.cbSize=0;
	 // Specify recording parameters
	 waveInOpen(&hWaveIn, WAVE_MAPPER,&pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
	 WAVEHDR      WaveInHdr;
	 // Set up and prepare header for input
	 WaveInHdr.lpData = (LPSTR)waveIn;
	 WaveInHdr.dwBufferLength = NUMPTS*2;
	 WaveInHdr.dwBytesRecorded=0;
	 WaveInHdr.dwUser = 0L;
	 WaveInHdr.dwFlags = 0L;
	 WaveInHdr.dwLoops = 0L;
	 waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	 HWAVEOUT hWaveOut;
	 cout << "playing..." << endl;
	 waveOutOpen(&hWaveOut, WAVE_MAPPER, &pFormat, 0, 0, WAVE_FORMAT_DIRECT);
	 waveOutWrite(hWaveOut, &WaveInHdr, sizeof(WaveInHdr)); // Playing the data
	 Sleep(3 * 1000); //Sleep for as long as there was recorded
	 waveInClose(hWaveIn);
	 waveOutClose(hWaveOut);
}

void StartRecord()
{
	 const int NUMPTS = 16025 * 3;   // 3 seconds
	 int sampleRate = 16025;  
	 // 'short int' is a 16-bit type; I request 16-bit samples below
	 // for 8-bit capture, you'd use 'unsigned char' or 'BYTE' 8-bit     types
	 HWAVEIN      hWaveIn;
	 MMRESULT result;
	 WAVEFORMATEX pFormat;
	 pFormat.wFormatTag=WAVE_FORMAT_PCM;     // simple, uncompressed format
	 pFormat.nChannels=1;                    //  1=mono, 2=stereo
	 pFormat.nSamplesPerSec=sampleRate;      // 8.0 kHz, 11.025 kHz, 22.05 kHz, and 44.1 kHz
	 pFormat.nAvgBytesPerSec=sampleRate*2;   // =  nSamplesPerSec × nBlockAlign
	 pFormat.nBlockAlign=2;                  // = (nChannels × wBitsPerSample) / 8
	 pFormat.wBitsPerSample=16;              //  16 for high quality, 8 for telephone-grade
	 pFormat.cbSize=0;
	 // Specify recording parameters
	 result = waveInOpen(&hWaveIn, WAVE_MAPPER,&pFormat, 0L, 0L, WAVE_FORMAT_DIRECT);
	 WAVEHDR      WaveInHdr;
	 // Set up and prepare header for input
	 WaveInHdr.lpData = (LPSTR)waveIn;
	 WaveInHdr.dwBufferLength = NUMPTS*2;
	 WaveInHdr.dwBytesRecorded=0;
	 WaveInHdr.dwUser = 0L;
	 WaveInHdr.dwFlags = 0L;
	 WaveInHdr.dwLoops = 0L;
	 waveInPrepareHeader(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	 // Insert a wave input buffer
	 result = waveInAddBuffer(hWaveIn, &WaveInHdr, sizeof(WAVEHDR));
	 // Commence sampling input
	 result = waveInStart(hWaveIn);
	 cout << "recording for 3 seconds..." << endl;
	 Sleep(1 * 1000);
	 cout<<"1s for initialization\n";
	 Sleep(2 * 1000);
	 // Wait until finished recording
	 waveInClose(hWaveIn);
	 PlayRecord();
}


int main()
{
	computeEnergyThresh();
	int e = fopen_s( &uni,"universe.txt", "w+");
	if(e != NULL)
	{
		printf("\nCannot open the file\n");
		system("pause");
		exit(1);
	}

	//training
	cout<<"----------------------------------------TRAINING--------------------------------------------\n";
	for(int i=0;i<10;i++){
		for(int j=1;j<=30;j++){
			char path[100];
			sprintf(path, "210101062_dataset/English/txt/210101062_E_%d_%d.txt", i, j);
			double x;
			int err;
			FILE *f1=NULL; // Read the data (x) from a text file
			cout<<path<<endl;
			err = fopen_s( &f1, (const char *)path, "r");
			if(err != NULL)
			{
				printf("\nCannot open the file\n");
				system("pause");
				exit(1);
			}                      
			num_samples=0;
			while( !feof(f1) )			//read the file
			{
				if( fscanf_s( f1, "%lf", &x) == 1)
				{
					samples[ num_samples ] = x;
					num_samples++;
				}
				else{
					//printf( "Error in line %d %d\n ", num_samples,file_no );
					char line[1000];
					fgets(line,100,f1);
				}
			}
			fclose(f1);
			getDCShift();
			getNFactor();
			for(int ns=0;ns<num_samples;ns++){
				samples[ns]=((samples[ns]-dcShift)*nFactor);
			}
			findSteadyFrames();
			int sampSize=0,curr=steadyStart,fr=0;
			while(curr<steadyEnd){
				int frame_size = (curr + framesize <= steadyEnd) ? framesize : steadyEnd - curr;
				for(int j=0; j<frame_size; j++){
					steady_frame[sampSize++] = samples[curr++];
				}
				computeRis(steady_frame+sampSize-frame_size,frame_size);
				computeAis();
				computeCis();
				for(int m=1;m<=P;m++)
					fprintf(uni,"%lf ",C[m]);
				fprintf(uni,"\n");
				fr++;
				//steady_frame[sampSize]*= 0.54-0.46*cos(2*PI*steady_frame[sampSize]/(framesize-1));
			}
			file_frames[file_no]=fr;
			file_no++;
			//cout<<"file "<<i<<" "<<j<<" completed\n";
		}
	}
	fclose(uni);

	ifstream file("universe.txt");
    string line, cell;
    if (!file.is_open()) {
        cerr << "Error: Could not open the file!" << endl;		    // Check if the file is open
        return 1;
    }
    while (getline(file, line)) {
        vector< double> row;				//vector to store one row of data
        stringstream lineStream(line);

        while (getline(lineStream, cell, ' ')) {
             double value;
            stringstream cellStream(cell);  
            cellStream >> value;			
            row.push_back(value);
        }
        universe.push_back(row);			// Add the row to the data vector
    }

	/*------------------------------Codebook generation can be commented once generated------------------------------------
	// Step 1: Start with a single vector representing the centroid of the universe
    codebook.push_back(calculateCentroid(universe));  // Initialize with the centroid of the universe

    // Set the convergence threshold (delta)
    double delta = 0.0001;

    // Step 2: Iteratively double the codebook size and refine using K-means
    while (codebook.size() < CB_SIZE) {
		cout<<"\n\ncode book of size "<<codebook.size()<<endl;
		cout<<"-----------------------------------------------------------------------------------------------------------------\n";
		for (int i=0;i< codebook.size();i++) {
			for (int j=0;j<codebook[i].size();j++) {
				cout << codebook[i][j] << " ";
			}
			cout << endl;
		}
		cout<<"-----------------------------------------------------------------------------------------------------------------\n\n";
        splitCodebook(codebook);						// Double the codebook size
        kMeans(universe, codebook, delta);	// Apply K-means to refine the centroids
    }

	FILE *fp = fopen("codebook.txt", "w");
    // Output the final codebook
    cout << "Final Codebook (by LBG) :" << endl;
	cout<<"---------------------------------------------------------------------------------------------------------------------\n";
    for (int i=0;i< codebook.size();i++) {
        for (int j=0;j<codebook[i].size();j++) {
            cout << codebook[i][j] << " ";
			fprintf(fp, "%lf  ", codebook[i][j]);
        }
		fprintf(fp,"\n");
        cout << endl;
    }
	cout<<"---------------------------------------------------------------------------------------------------------------------\n";
	fclose(fp);*/
    // Close the file
    file.close();


	ifstream file1("codebook.txt"); // Open the CSV file
    if (!file1.is_open()) {
        cerr << "Error: Could not open the file!" << endl;
        return 1;
    }
    while (getline(file1, line)) {
        vector< double> row;				// Vector to store one row of data
        stringstream lineStream(line);

        while (getline(lineStream, cell, ' ')) {
             double value;
            stringstream cellStream(cell);  
            cellStream >> value;			
            row.push_back(value);			
        }
        codebook.push_back(row);			// Add the row to the data vector
    }
	file_no=0;
	int frames=0;

	for(int i=0;i<10;i++){
		for(int j=1;j<=TRAIN_FILES;j++){
			resetModel();
			T=file_frames[file_no];
			for(int k=1;k<=T;k++){
				double temp_ci[13]={0};
				for(int r=1;r<=12;r++) temp_ci[r]=universe[frames][r-1];
				O[k] = genObservations(temp_ci);
				frames++;
			}

			int iteration=0;
			p_star = 0, prev_p_star = -1;
			while(p_star > prev_p_star && iteration < 1000){
				iteration++;
				prev_p_star = p_star; 
				forwardProcedure();
				backwardProcedure();
				viterbiAlgorithm();
				computeXi();
				computeGamma();
				reestimateParameters();
			}
			addModels();
			file_no++;
		}
		computeAvgModel(TRAIN_FILES);
		writeAvgModels(i);
		resetAvgModel();
	}

	cout<<"MENU:\n1. Test recordings\n2. Live testing\nEnter choice(1/2): ";
	int choice;
	cin>>choice;
	if(choice == 1){
		testRecordings();
	} 
	else if(choice == 2){
		StartRecord();
		testingLive();
	} 
	else{
		cout<<"Invalid choice!\n";
	}
	return 0;
}