#include "abc.h"

/* constant values */
using namespace std;
using namespace std::chrono;
const double PI_VALUE = 3.14159265358979323846;
constexpr double inf = std::numeric_limits<double>::infinity();
const double E_VALUE = 2.7182818284590452354;
int iterations = 0;


/* bee structure */
struct Bee {
	vector<int> pos;
	double cost = 0.0f;
};

int N = 32;
int M = 32;

void multiply_matrix(vector<vector<double>> &A, vector<double> &B, vector<double> &C) {
	// A - matrix, B - array, C - result
	for (int j = 0; j < N; ++j)
		C[j] = 0;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			C[i] += A[j][i] * B[j];
		}
	}
}

vector<vector<double>> K_t, K;

void transpose(vector<vector<double>> &A) {
	double temp;
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < i; ++j) {
			temp = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = temp;
		}
	}
}


/* class that holds the optimizable functions ready for use */
class OptimizableFunction {
public:
	/* bounds for the used function */
	int m_lowerBound;
	int m_upperBound;
	vector<double> F;
	vector<double> T;

	/* constructor */
	OptimizableFunction(vector<double> &values) {
		setBounds();
		F = values;
		T.resize(N);
		multiply_matrix(K_t, F, T);
	}

	/* set the bounds depending on the function */
	void setBounds() {
        // границы для значений координат
        m_lowerBound = 0;
        m_upperBound = N-1;
	}

	/* return function value for a vector of positions */
	double getResult(vector<int> &numbers) {
        // here the error after using the three chosen coefficients must be calculated - compress and decompress
        // set to zero all coefficients except 3 picked by the bee, perform reverse Tchebichef transform and calculate MSE
        iterations--;
		for (int x: numbers) {
			if (x < 0 || x > N-1)
				return inf;
		}
		vector<double> T_reserved(N, 0.);
        for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M; ++j) {
				if (numbers[j] == i) {
					T_reserved[i] = T[i];
					break;
				}
			}
		}
		vector<double> R(N, 0.);
		multiply_matrix(K, T_reserved, R);
		double sum = 0.;
		for (int i = 0; i < N; ++i) {
			sum += (F[i] - R[i]) * (F[i] - R[i]);
		}
		sum /= N;
		return sum;
	}
};

/* returns the cumulative sum of a vector of numbers */
vector<double> cumSum(vector<double> P) {
	double partialSum = 0;
	for (int i = 0; i < P.size(); i++) {
		partialSum += P[i];
		P[i] = partialSum;
	}
	return P;
}

/* returns index of potentially useful solution */
int fitnessProportionateSelection(vector<double> P) {
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> randDist(0, 1);
	double randNumber = randDist(gen);
	vector<double> cumVector = cumSum(P);
	for (int i = 0; i < cumVector.size(); i++) {
		if (randNumber <= cumVector[i]) {
			return i;
		}
	}
	return NULL;
}

double* solve(double* numbers, double* chebyshev, int N_, int M_) {
	cout << "kek" << endl;
	N = N_;
	M = M_;
	K.resize(N, vector<double>(N));
	K_t.resize(N, vector<double>(N));
	for (int i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			K[i][j] = chebyshev[i * N + j];
			K_t[j][i] = chebyshev[i * N + j];
		}
	}
	//transpose(K_t);
	/* define infinity */
	/* define random generator template */
	random_device rd;
	mt19937 gen(rd());
	/* define old random generator with current time as seed */
	srand((unsigned int)time(NULL));

	/* get through all functions */
	/*get through all distance sizes == размер вектора pos (в моем случае 32)*/
	int distanceSize = M;
    auto start = high_resolution_clock::now();

	/* maximum number of experiments */
	int nrOfExperiments = 1;
	/* maximum number of iterations */
	int maxIterations = 6000;

	/* vector for holding the best cost of each experiment*/
	vector<double> bestCostExperiments(nrOfExperiments, 0);
	/* mean cost of all experiments */
	double allMeanCost = 0.0f;
	/* standard deviation cost of all experiments */
	double allSDCost = 0.0f;

	/* create a function object to be tested; parameter is used for choosing a certain function */
	vector<double> otch;
	otch.assign(numbers, numbers+N);
	//vector<double> otch = {995, 995, 995, 995, 995, 995, 995, 995, 1000, 997, 995, 994, 992, 993, 992, 989, 988, 987, 990, 993, 989, 988, 986, 988, 993, 997, 993, 986, 983, 977, 979, 975};
	OptimizableFunction optimizableFunction(otch);
	/* get the bound values for the function */
	int distanceMin = optimizableFunction.m_lowerBound;
	int distanceMax = optimizableFunction.m_upperBound;

	/* colony size of bees (or population size) */
	int beesPopulation = 125;
	/* employed bees (sending them onto the food sources to measure their nectar amounts / fitness value) */
	int beesEmployed = int(50.0f / 100.0f * beesPopulation);
	/* onlooker bees (select the food sources using the nectar information / fitness value) */
	int beesOnlooker = beesEmployed;
	/* scout bees (sent to the selected food sources) */
	int beesScout = 1;
	/* limit used for determining the worth of finding a food source */
	int limit = int(round(0.6 * distanceSize * beesPopulation));
	/* acceleration parameter used for finding potential food sources faster */
	int accel = 1;

	/* loop through all experiments */
	for (int experiment = 0; experiment < nrOfExperiments; experiment++) {
		/* vector for the entire population of bees */
		vector<Bee> beesVector;
		/* initialize the vector with empty bees */
		for (int i = 0; i < beesPopulation; i++) {
			Bee initBee;
			beesVector.push_back(initBee);
		}
		/* initialize the best solution with the worst cost */
		Bee bestSolution;
		bestSolution.cost = inf;
		/* initializing distances for the bees vector with random values and calculating costs */
		for (int i = 0; i < beesPopulation; i++) {
			vector<int> indeces(N);
			for (int j = 0; j < N; ++j)
				indeces[j] = j;
			for (int j = 0; j < distanceSize; j++) {
				int chosen = rand() % indeces.size();
				beesVector[i].pos.push_back(indeces[chosen]);
				indeces.erase(indeces.begin() + chosen);
			}
			beesVector[i].cost = optimizableFunction.getResult(beesVector[i].pos);
			if (beesVector[i].cost <= bestSolution.cost) {
				bestSolution = beesVector[i];
			}
		}
		cout << beesVector[0].cost << endl;
		/* vector for counting the abandoned bees */
		vector<double> abandonedBees(beesPopulation, 0);
		/* vector for keeping in memory the best cost for every iteration */
		vector<double> bestCost(maxIterations, 0);
        iterations = maxIterations;
		/* loop through all iterations */
        int it = 0;
		while (iterations > 0) {

			/* employeed bees phase */
			for (int i = 0; i < beesEmployed; i++) {

				/* choose a random bee that isn't the current one */
				vector<int> randomBees;
				for (int j = 0; j < beesEmployed; j++) {
					if (j != i) {
						randomBees.push_back(j);
					}
				}
				int randomBeesIndex = randomBees[rand() % randomBees.size()];

				/* calculate a different acceleration coefficient for every distance */
				vector<double> accelCoef;
				uniform_real_distribution<> acc(-1, +1);
				for (int j = 0; j < distanceSize; j++) {
					accelCoef.push_back(accel * acc(gen));
				}

				/* define a new bee */
				Bee newBee;
				/* get the new bee position */
				for (int j = 0; j < distanceSize; j++) {
					/* new bee position is equal to current bee's position + (current bee's position - random bee's position) * acceleration coefficient */
					newBee.pos.push_back(max(0., min((double)distanceMax, beesVector[i].pos[j] + (beesVector[i].pos[j] - beesVector[randomBeesIndex].pos[j]) * accelCoef[j])));
				}
				/* calculate the new cost */
				newBee.cost = optimizableFunction.getResult(newBee.pos);

				/* if the new cost is better */
				if (newBee.cost <= beesVector[i].cost) {
					/* replace the old cost with the better one */
					beesVector[i] = newBee;
					abandonedBees[i] = 0;
				}
				/* else abandon it */
				else {
					abandonedBees[i] += 1;
				}
			}

			/* vector for calculating the fitness values */
			vector<double> fitnessValues(beesEmployed, 0);
			/* sum of fitness values */
			double fSum = 0;
			/* average cost of all the bees */
			double averageCost = 0;
			/* calculating the average cost */
			for (int i = 0; i < fitnessValues.size(); i++) {
				averageCost += beesVector[i].cost;
			}
			averageCost /= fitnessValues.size();
			/* calculating the fitness values */
			for (int i = 0; i < fitnessValues.size(); i++) {
				fitnessValues[i] = pow(E_VALUE, -beesVector[i].cost / averageCost);
				fSum += fitnessValues[i];
			}
			/* calculating probability of being selected */
			vector<double> prob(beesEmployed, 0);
			for (int i = 0; i < prob.size(); i++) {
				prob[i] = fitnessValues[i] / fSum;
			}

			/* onlooker bees phase */
			for (int m = 0; m < beesOnlooker; m++) {

				/* select a bee index based on the probability and make it the current one */
				int i = fitnessProportionateSelection(prob);

				/* choose a random bee that isn't the current one */
				vector<int> randomBees;
				for (int j = 0; j < beesEmployed; j++) {
					if (j != i) {
						randomBees.push_back(j);
					}
				}
				int randomBeesIndex = randomBees[rand() % randomBees.size()];

				/* calculate a different acceleration coefficient for every distance */
				vector<double> accelCoef;
				uniform_real_distribution<> acc(-1, +1);
				for (int j = 0; j < distanceSize; j++) {
					accelCoef.push_back(accel * acc(gen));
				}

				/* define a new bee */
				Bee newBee;
				/* get the new bee position */
				for (int j = 0; j < distanceSize; j++) {
					/* new bee position is equal to current bee's position + (current bee's position - random bee's position) * acceleration coefficient */
					newBee.pos.push_back(max(0., min((double)distanceMax, beesVector[i].pos[j] + (beesVector[i].pos[j] - beesVector[randomBeesIndex].pos[j]) * accelCoef[j])));
				}
				/* calculate the new cost */
				newBee.cost = optimizableFunction.getResult(newBee.pos);

				/* if the new cost is better */
				if (newBee.cost <= beesVector[i].cost) {
					/* replace the old cost with the better one */
					beesVector[i] = newBee;
					abandonedBees[i] = 0;
				}
				/* else abandon it */
				else {
					abandonedBees[i] += 1;
				}
			}

			/* scout bees phase */
			for (int itScout = 0; itScout < beesScout; itScout++) {
				for (int i = 0; i < beesEmployed; i++) {
					/* if the abandoned bee is over the abandonment limit, make it search again */
					if (abandonedBees[i] >= limit) {
						uniform_real_distribution<> distance(distanceMin, distanceMax);
						vector<int> indeces(N);
						for (int j = 0; j < N; ++j)
							indeces[j] = j;
						for (int j = 0; j < distanceSize; j++) {
							int chosen = rand() % indeces.size();
							beesVector[i].pos.push_back(indeces[chosen]);
							indeces.erase(indeces.begin() + chosen);
						}
						beesVector[i].cost = optimizableFunction.getResult(beesVector[i].pos);
						abandonedBees[i] = 0;
					}
				}
			}

			/* calculating the best cost */
			for (int i = 0; i < beesEmployed; i++) {
				if (beesVector[i].cost <= bestSolution.cost) {
					bestSolution = beesVector[i];
				}
			}
			/* saving the best cost in a vector */
			bestCost[it] = bestSolution.cost;

			if (iterations <= 0) {
				bestCostExperiments[experiment] = bestCost[it];
				//cout << "Best cost for experiment " << experiment << ": " << bestCostExperiments[experiment] << endl;
				//int* res_arr = new int[M];
				//for (int l = 0; l < M; ++l)
				//	res_arr[l] = bestSolution.pos[l];
				double* res_arr = new double[N];
				vector<double> F = otch;
				vector<double> T(N);
				multiply_matrix(K_t, F, T);
				vector<double> T_reserved(N, 0.);
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < M; ++j) {
						if (bestSolution.pos[j] == i) {
							T_reserved[i] = T[i];
							break;
						}
					}
				}
				vector<double> R(N, 0.);
				multiply_matrix(K, T_reserved, R);
				copy(R.begin(), R.end(), res_arr);
				return res_arr;
			}

		}
	}
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    //cout << duration.count() << endl;
}