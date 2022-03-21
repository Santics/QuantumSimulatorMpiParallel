#include"Dependence.h"
#include"Gate.h"
#include"State.h"

#define QUBITNUM 20

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	int numprocs, myid;
	int nonlocalQubitNum, localQubitNum;
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	double probability[2] = { 0 };
	double* multiProbability = NULL;
	double specificProbability = 0;

	nonlocalQubitNum = log(numprocs) / log(2);
	localQubitNum = QUBITNUM - nonlocalQubitNum;

	//Gate Intialization
	Matrix2cd mI, mX, mY, mZ, mH, mS, mT;

	mI << 1, 0,
		0, 1;

	mX << 0, 1,
		1, 0;

	mY << 0, complex<double> {0, -1},
		complex<double>{ 0, 1 }, 0;

	mZ << 1, 0,
		0, -1;

	mH << 1 / sqrt(2), 1 / sqrt(2),
		1 / sqrt(2), -1 / sqrt(2);

	mS << 1, 0,
		0, complex<double>{ 0, 1 };

	mT << 1, 0,
		0, exp(complex<double>{ 0, M_PI / 4 });

	//default gate intialize
	Gate I(mI, 1), X(mX, 1), Y(mY, 1), Z(mZ, 1), H(mH, 1), S(mS, 1), T(mT, 1);
	vector<Gate> GateSet;
	GateSet.push_back(I);
	GateSet.push_back(X);
	GateSet.push_back(Y);
	GateSet.push_back(Z);
	GateSet.push_back(H);
	GateSet.push_back(S);
	GateSet.push_back(T);


	//generalQubitGate Test part
	//MatrixXcd mDemo = MatrixXcd::Random(4,4);
	Matrix4cd mDemo;

	mDemo << 1, 0, 0, 0
		, 0, 1, 0, 0
		, 0, 0, 1 / sqrt(2), 1 / sqrt(2)
		, 0, 0, 1 / sqrt(2), -1 / sqrt(2);

	Gate Demo(mDemo, 2);
	vector<int> qubitSequence = { 4,3 };
	vector<int> measureSequence = { 2,3,4 };


	//edit this part manually to implement some specific quantum algorithms
	//just give an array of datatype:Gate(see "Gate.h" for more infomation)
	//GUI will be added in the following version to edit the circuit directly

	State a(QUBITNUM);
	a.Initialize("11111111111111111111");

	a.Permutation(5, 1);
	//a.Permutation(4, 2);
	//a.Permutation(0, 1);
	//a.Permutation(3, 0);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 3, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 4);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 4, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	Demo(localQubitNum, &a, qubitSequence);
	Demo(localQubitNum, &a, qubitSequence);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 3, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 4);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 4, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	Demo(localQubitNum, &a, qubitSequence);
	Demo(localQubitNum, &a, qubitSequence);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 3, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 4);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 4, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	Demo(localQubitNum, &a, qubitSequence);
	Demo(localQubitNum, &a, qubitSequence);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 3, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	GateSet[4](localQubitNum, &a, 0);
	GateSet[4](localQubitNum, &a, 5);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 4);
	GateSet[4](localQubitNum, &a, 3);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 2);
	GateSet[4](localQubitNum, &a, 1);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 4);
	GateSet[4](localQubitNum, &a, 1, 4, 2);
	GateSet[4](localQubitNum, &a, 1, 3, 5);
	GateSet[4](localQubitNum, &a, 1, 5, 0);
	Demo(localQubitNum, &a, qubitSequence);
	Demo(localQubitNum, &a, qubitSequence);

	//a.measure(localQubitNum,3,probability);
	multiProbability = a.measure(localQubitNum, measureSequence);


	if (myid == 0)
	{
		cout << "Probability 000 :" << multiProbability[0] << endl;
		cout << "Probability 001 :" << multiProbability[1] << endl;
		cout << "Probability 010 :" << multiProbability[2] << endl;
		cout << "Probability 011 :" << multiProbability[3] << endl;
		cout << "Probability 100 :" << multiProbability[4] << endl;
		cout << "Probability 101 :" << multiProbability[5] << endl;
		cout << "Probability 110 :" << multiProbability[6] << endl;
		cout << "Probability 111 :" << multiProbability[7] << endl;
	}

	//if (myid == 0)
	//{
	//	cout << "Probability 00 :" << multiProbability[0] << endl;
	//	cout << "Probability 01 :" << multiProbability[1] << endl;
	//	cout << "Probability 10 :" << multiProbability[2] << endl;
	//	cout << "Probability 11 :" << multiProbability[3] << endl;
	//}

	//if (myid == 0)
	//{
	//	cout << "Probability 0 :" << multiProbability[0] << endl;
	//	cout << "Probability 1 :" << multiProbability[1] << endl;
	//}

	a.measure(localQubitNum, measureSequence, "000", specificProbability);
	if (myid == 0)
		cout << endl << "specific: " << specificProbability << endl;

	//if (myid == 0)
	//{
	//	cout << "Probability 0 :" << probability[0] << endl;
	//	cout << "Probability 1 :" << probability[1] << endl;
	//}

	//when the algorithm is done, call python to print the circuit
	MPI_Finalize();
}