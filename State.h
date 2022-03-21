#pragma once
#include"Dependence.h"

class State
{
	/****************/
	//Attribute
	/****************/
public:

	//the number of qubits
	int qubitNum;
	//array denoting the sequence of local qubits
	int* nonLocalQubits;
	int* localQubits;
	//data of the state vector
	VectorXcd stateVector;


	/****************/
	//Member Function
	/****************/


//developer based functions
private:
	//judge whether a specific qubit is local
	bool whetherQubitLocal(int localNum, int qubit, int& orderofAimQubit);
	bool whetherQubitNonlocal(int nonlocalNum, int qubit, int& orderofAimQubit);
	void LocateQubit(int nonlocalQubit, int localQubit, int nonlocalQubitNum, int localQubitNum, int& nonlocalIndex, int& localIndex);
	int qubitFind(vector<int> qubits, int  targetQubit);

	//user based functions
public:
	friend class Gate;
	//constructor function
	State(int localQbitNum);

	void Initialize(std::string initialQubitString);

	void Permutation(int nonlocalQubit, int localQubit);
	void Permutation(int permutationNum, int* nonlocalQubit, int* localQubit);

	void Permutation(int nonlocalQubit1, int nonlocalQubit2, int localQubit1, int localQubit2);


	void measure(int localNum, int measureQubit, double* probability);
	double* measure(int localNum, vector<int>& measureQubit);
	void measure(int localNum, vector<int>& measureQubit, std::string qubits, double& probability);
};

