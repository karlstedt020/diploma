#include "abc.h"

/* constant values */
using namespace std;
using namespace std::chrono;
const double PI_VALUE = 3.14159265358979323846;
constexpr double inf = std::numeric_limits<double>::infinity();
const double E_VALUE = 2.7182818284590452354;
int iterations = 0;
double chebyshevT(int n, double x) {
    if(n == 0)
        return 1.0;
    if(n == 1)
        return x;
    double Tnm2 = 1.0;  // T0(x)
    double Tnm1 = x;    // T1(x)
    double Tn_val = 0.0;
    for (int k = 2; k <= n; k++) {
        Tn_val = 2.0 * x * Tnm1 - Tnm2;
        Tnm2 = Tnm1;
        Tnm1 = Tn_val;
    }
    return Tn_val;
}

// Build the Chebyshev Vandermonde matrix A such that A[j][n] = T_n(x_j)
// where x_j are the sample nodes.
vector<vector<double>> computeChebyshevMatrix(const vector<double>& x) {
    int m = x.size();
    vector<vector<double>> A(m, vector<double>(m, 0.0));
    for (int j = 0; j < m; j++) {
        for (int n = 0; n < m; n++) {
            A[j][n] = chebyshevT(n, x[j]);
        }
    }
    return A;
}

// Solve the linear system A * c = b using Gaussian elimination.
// (A is assumed to be square and invertible.)
vector<double> solveLinearSystem(vector<vector<double>> A, vector<double> b) {
    int m = A.size();
    
    // Forward elimination with partial pivoting.
    for (int k = 0; k < m; k++) {
        // Find pivot row
        int pivot = k;
        double maxVal = fabs(A[k][k]);
        for (int i = k+1; i < m; i++) {
            if (fabs(A[i][k]) > maxVal) {
                maxVal = fabs(A[i][k]);
                pivot = i;
            }
        }
        if (pivot != k) {
            swap(A[k], A[pivot]);
            swap(b[k], b[pivot]);
        }
        // Eliminate entries below the pivot.
        for (int i = k+1; i < m; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < m; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    
    // Back substitution.
    vector<double> c(m, 0.0);
    for (int i = m - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < m; j++) {
            sum += A[i][j] * c[j];
        }
        c[i] = (b[i] - sum) / A[i][i];
    }
    return c;
}

/* bee structure */
struct Bee {
	vector<int> pos;
	double cost = 0.0f;
};

const int N = 32;
const int M = 4;

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
	vector<double> x;

	/* constructor */
	OptimizableFunction(vector<double> &values) {
		setBounds();
		F = values;
		T.resize(N);
		x.resize(N, 0.0);
		for (int j = 0; j < N; j++) {
			x[j] = cos(M_PI * j / (N - 1));
		}
		vector<vector<double>> A = computeChebyshevMatrix(x);
		T = solveLinearSystem(A, values);
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
		vector<double> T_reserved(N, 0.0);
		for (int j = 0; j < N; ++j) {
			for (int i = 0; i < M; ++i) {
				if (j == numbers[i]) {
					T_reserved[j] = T[j];
					break;
				}
			}
		}
		vector<double> ecg_rec(N, 0.0);
		for (int j = 0; j < N; j++) {
			// Option 1: Direct summation using recurrence.
			double sum = 0.0;
			for (int n = 0; n < N; n++) {
				sum += T_reserved[n] * chebyshevT(n, x[j]);
			}
			ecg_rec[j] = sum;
		}
		double sum = 0.;
		for (int i = 0; i < N; ++i) {
			sum += (F[i] - ecg_rec[i]) * (F[i] - ecg_rec[i]);
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
	//int limit = int(round(0.6 * distanceSize * beesPopulation));
	int limit = 20;
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
				vector<double> ecg_rec(N, 0.0);
				vector<double> x;
				vector<double> T_reserved(N, 0.0);
				x.resize(N, 0.0);
				for (int j = 0; j < N; j++) {
					x[j] = cos(M_PI * j / (N - 1));
				}
				vector<vector<double>> A = computeChebyshevMatrix(x);
				T = solveLinearSystem(A, F);
				for (int j = 0; j < N; ++j) {
					for (int i = 0; i < M; ++i) {
						if (j == bestSolution.pos[i]) {
							T_reserved[j] = T[j];
							break;
						}
					}
				}
				for (int j = 0; j < N; j++) {
					// Option 1: Direct summation using recurrence.
					double sum = 0.0;
					for (int n = 0; n < N; n++) {
						sum += T_reserved[n] * chebyshevT(n, x[j]);
					}
					ecg_rec[j] = sum;
				}
				copy(ecg_rec.begin(), ecg_rec.end(), res_arr);
				double sum = 0.;
				for (int i = 0; i < N; ++i) {
					//sum += (F[i] - ecg_rec[i]) * (F[i] - ecg_rec[i]);
					sum += std::fabs((F[i] - ecg_rec[i]) / F[i]);
				}
				sum /= N;
				cout << sum << endl;
				return res_arr;
			}

		}
	}
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    //cout << duration.count() << endl;
}