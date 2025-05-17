AIM: Digit Recognition
EXECUTION: Build and run the code on Visual Studio 2010. Use F5 key to run the code.

INPUT: 210101062_dataset
OUTPUT: Digit prediction on a new recording.

CONSTANTS:
1. M: The number of observation symbols per state (32)
2. N: The number of states (5)
3. MAX_T: Maximum number of observations (160)
4. DELTA: Used to terminate the K-means algorithm (0.0001)
5. EPSILON: The splitting parameter (0.03)
6. tokhuraWeights: Array to store the 12 Tokhura weight values.
7. P: The number of past samples to consider (12)
8. CB_SIZE: Size of the codebook (32)

VARIABLES:
1. A[N+1][N+1]: Array to store the state transition probability distribution.
2. B[N+1][M+1]: Array to store the observation symbol probability distribution.
3. Pi[N+1]: Array to store the initial state distribution.
4. T: Variable to store number of observations
5. O[MAX_T+1]: Array to store the observation sequence.

6. alpha[MAX_T+1][N+1]: Array to store the forward variable values.
7. beta[MAX_T+1][N+1]: Array to store the backward variable values.
8. P_O_given_lambda: Long double variable to store the solution of problem 1.

9. delta[MAX_T+1][N+1]: Array to store the best score along a single path, at time t, which accounts for the first t observations and ends in state i.
10. psi[MAX_T+1][N+1]: Array to keep track of the argument that maximized delta[t+1][j] for each t and j.
11. q_star[MAX_T+1]: Array to store the optimal state sequence.
12. p_star: long double to store the solution of problem-2.

13. gamma[MAX_T+1][N+1]: Array to store the gamma values.
14. xi[MAX_T+1][N+1][N+1]: Array to store the probability of being in state Si at time t & Sj at time t+1 given the observation sequence & the model.
15. A_bar[N+1][N+1]: Array to store the reestimated state transition probability distribution.
16. B_bar[N+1][M+1]: Array to store the reestimated observation symbol probability distribution.
17. Pi_bar[N+1]: Array to store the reestimated initial state distribution.


FUNCTIONS:
1. tokhuraDistance():
	Function to calculate the Tokhura distance between the current centroid and current row of universe.
2. kMeans():
	Function to implement Lloyd's algorithm for vector quantization.
3. splitCodebook():
	Function to double the size of the codebook by splitting each current code vector.
4. forwardProcedure():
	Function to inductively solve for the forward variable, where alpha is the probability of the partial observation sequence until time t and state Si at time t given model lambda. 
5. backwardProcedure():
	Function to inductively solve for the backward variable, where beta is the probability of the partial observation from t+1 till the end, given state Si at time t and model lambda. 
6. viterbiAlgorithm():
	Function to find the single best state sequence for the given observation sequence.
7. computeXi():
	Function to compute the probability of being in state Si at time t & Sj at time t+1 given the observation sequence & the model.
8. computeGamma():
	Function to compute the probability of being in state Si at time t given the observation sequence & the model.
9. reestimateParameters():
	Function to iteratively update and improve the HMM parameters.
10. genObservations():
    Function to generate observation sequence from the ci's.
11. addModels():
    Function to add models to find average.
12. computeAvgModel():
    Function to compute average model.
13. testingLive():
    Function for testing a live recording.
14. testRecordings():
    Function for testing the test recordings.
15. resetModel():
    Function to reset the current model to the unbiased/ergodic model.
16. resetAvgModel():
    Function to reset the average model to zeroes.
17. computeEnergyThresh():
	Function to calculate the silence energy.
18. findSteadyFrames():
	Function to find the steady frames.
19. computeRis():
	Function to compute the Ri values.
20. computeAis():
	Function to compute the Ai values.
21. computeCis():
	Function to compute the Ci values.

PROCEDURE:
1. Generate Ci vector for each recording file and dump into universe file.
2. Use this universe file to create the codebook using the LBG code.
3. Using this codebook generate observation sequence for each file.
4. Obtain observation seq, A Matrix (5x5) , B Matrix (5x32) and Pi Matrix (1x5)
5. Use HMM Problem 1 , 2 and 3 to improve the model i.e. A matrix and B matrix
