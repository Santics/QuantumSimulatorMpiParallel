#pragma once
#include"Dependence.h"
#include"State.h"

class Gate
{
	/****************/
	//Attribute
	/****************/
private:

	int dim;
	MatrixXcd gateMatrix;


	/****************/
	//Member Function
	/****************/

private:
	//auxiliary functions
	//developer based, not used based

	bool searchElement(int target, int arrayLength, vector<int>& array);
	void indiceBinaryGenerator(vector<int>& indice, vector<int>& targetSequence, int dim);

	//
	void singleQubitGateOperates(int localNum, int orderofAimQubit, State* state);
	void controlledGateOperates(int localNum, bool controlGateType, int orderofControlQubit, int orderofTargetQubit, State* state);
	void generalQubitGateOperates(int localNum, vector<int>& qubitSequence, State* state);



public:

	//constructor function
	//when gateType = 0 , size of attribute gateMatrix is 2*2 :general single qubit gate or controlled gate;
	//when gateType = 1 , size of attribute gateMatrix is 4*4 :general two qubits gate;
	Gate(MatrixXcd gateMatrix, int qubitSize);



	//overloaded operator "()" to let Gate operates
	//you can ues syntax like this
	//someSingleQbitGate(qubit)  or  someCGate(controlQubit,targetQubit)
	//to let some specific gates operates on some specific qubits
	void operator()(int localNum, State* state, int singleGateQubit);
	//control gate
	//controlQubitType indicates the type of control qubit
	//0-control or 1-control; 
	void operator()(int localNum, State* state, bool controlQubitType, int controlQubit, int targetQubit);
	//general qubits gate
	void operator()(int localNum, State* state, vector<int>& qubitSequence);
};

