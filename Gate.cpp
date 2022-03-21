#include "Gate.h"

//*****************************************************
//developer based auxiliary functions

//judge whether a specific qubit is local
bool Gate::searchElement(int target, int arrayLength, vector<int>& array)
{
	int i = 0;
	bool flag = false;
	for (i = 0; i < arrayLength; i++)
	{
		if (target == array[i])
		{
			flag = true;
			break;
		}
	}

	return flag;
}

void Gate::indiceBinaryGenerator(vector<int>& indice, vector<int>& targetSequence, int dim)
{
	int i = 0, j = 0, k = 0;


	for (i = 0; i < pow(2, dim); i++)
	{
		indice[i] = 0;
		for (j = 1, k = 0; k < dim; j = j << 1, k++)
		{
			if (i & j)
			{
				indice[i] += pow(2, targetSequence[k]);
			}
		}

	}

	return;
}

void Gate::singleQubitGateOperates(int localNum, int orderofAimQubit, State* state)
{
	long long int i;
	long long int gap = pow(2, orderofAimQubit);
	std::complex<double> tmp1, tmp2;

	for (i = 0; i < pow(2, localNum); i++)
	{
		tmp1 = state->stateVector[i];
		tmp2 = state->stateVector[i + gap];

		state->stateVector[i] = tmp1 * this->gateMatrix(0, 0) + tmp2 * this->gateMatrix(0, 1);
		state->stateVector[i + gap] = tmp1 * this->gateMatrix(1, 0) + tmp2 * this->gateMatrix(1, 1);

		if (((i + 1) % gap) == 0)
			i += gap;
	}

	return;
}

void Gate::controlledGateOperates(int localNum, bool controlGateType, int orderofControlQubit, int orderofTargetQubit, State* state)
{
	long long int i = 0, j = 0, baseIndice = 0;
	int k = 0;
	vector<int> tempArray = { orderofTargetQubit,orderofControlQubit };
	vector<int> res(localNum - 2), indice(4);
	Vector2cd temp;

	for (i = 0, j = 0; i < localNum; i++)
	{
		if (!searchElement(i, 2, tempArray))
		{
			res[j] = i;
			j++;
		}
	}
	indiceBinaryGenerator(indice, tempArray, 2);
	for (i = 0; i < pow(2, localNum - 2); i++)
	{
		baseIndice = 0;
		for (j = 1, k = 0; k < localNum - 2; j = j << 1, k++)
		{
			if (i & j)
				baseIndice += pow(2, res[k]);
		}
		if (controlGateType)
		{
			temp[0] = state->stateVector[baseIndice + indice[2]];
			temp[1] = state->stateVector[baseIndice + indice[3]];
		}
		else
		{
			temp[0] = state->stateVector[baseIndice + indice[0]];
			temp[1] = state->stateVector[baseIndice + indice[1]];
		}

		temp = this->gateMatrix * temp;

		if (controlGateType)
		{
			state->stateVector[baseIndice + indice[2]] = temp[0];
			state->stateVector[baseIndice + indice[3]] = temp[1];
		}
		else
		{
			state->stateVector[baseIndice + indice[0]] = temp[0];
			state->stateVector[baseIndice + indice[1]] = temp[1];
		}
	}

	return;
}

void Gate::generalQubitGateOperates(int localNum, vector<int>& qubitOrderSequence, State* state)
{
	long long int i = 0, k = 0, baseIndice = 0, indice = 0;
	int j = 0, m = 0, n = 0;
	vector<int>res(localNum - this->dim);
	VectorXcd tmp = VectorXcd::Zero(pow(2, this->dim));


	for (i = 0, j = 0; i < localNum; i++)
	{
		if (!searchElement(i, this->dim, qubitOrderSequence))
		{
			res[j] = i;
			j++;

		}
	}

	for (i = 0; i < pow(2, localNum - this->dim); i++)
	{
		baseIndice = 0;
		for (j = 1, m = 0; m < (localNum - this->dim); j = j << 1, m++)
		{
			if (i & j)
			{
				baseIndice += pow(2, res[m]);
			}
		}
		for (k = 0; k < pow(2, this->dim); k++)
		{
			indice = baseIndice;
			for (m = 1, n = 0; n < this->dim; m = m << 1, n++)
			{
				if (k & m)
				{
					indice += pow(2, qubitOrderSequence[n]);
				}
			}
			tmp[k] = state->stateVector[indice];
		}

		tmp = this->gateMatrix * tmp;

		for (k = 0; k < pow(2, this->dim); k++)
		{
			indice = baseIndice;
			for (m = 1, n = 0; n < this->dim; m = m << 1, n++)
			{
				if (k & m)
				{
					indice += pow(2, qubitOrderSequence[n]);
				}
			}
			state->stateVector[indice] = tmp[k];
		}
	}

	return;
}



//*****************************************
//user based functions
Gate::Gate(MatrixXcd gateMatrix, int qubitSize)
{
	int	i = 0, j = 0;
	long long int gateSize = pow(2, qubitSize);
	this->dim = qubitSize;
	this->gateMatrix = MatrixXcd::Zero(gateSize, gateSize);

	for (i = 0; i < gateSize; i++)
	{
		for (j = 0; j < gateSize; j++)
		{
			this->gateMatrix(i, j) = gateMatrix(i, j);
		}
	}
	return;
}


void Gate::operator()(int localNum, State* state, int singleGateQubit)
{
	bool qubitLocal = 0;
	int orderofAimQubit;
	int partitionNum, partitionElementNum;

	//if qubit to be operated is nonlocal
	//just make it local!
	qubitLocal = state->whetherQubitLocal(localNum, singleGateQubit, orderofAimQubit);
	if (!qubitLocal)
	{
		state->Permutation(singleGateQubit, state->localQubits[localNum - 1]);
		orderofAimQubit = localNum - 1;
	}

	//Gate operates on single qubit
	singleQubitGateOperates(localNum, orderofAimQubit, state);


	return;
}

void Gate::operator()(int localNum, State* state, bool controlQubitType, int controlQubit, int targetQubit)
{
	bool controlQubitLocal, targetQubitLocal;
	int orderofControlQubit, orderofTargetQubit;
	int orderofNonlocalControlQubit;
	int myid = 0;

	controlQubitLocal = state->whetherQubitLocal(localNum, controlQubit, orderofControlQubit);
	targetQubitLocal = state->whetherQubitLocal(localNum, targetQubit, orderofTargetQubit);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (!controlQubitLocal && !targetQubitLocal)
	{
		state->Permutation(targetQubit, state->localQubits[localNum - 1]);
		orderofTargetQubit = localNum - 1;
		state->whetherQubitNonlocal(state->qubitNum - localNum, controlQubit, orderofNonlocalControlQubit);
		if (!(((myid >> orderofNonlocalControlQubit) % 2) ^ (controlQubitType)))
			singleQubitGateOperates(localNum, orderofTargetQubit, state);
		return;
	}
	else if (controlQubitLocal && !targetQubitLocal)
	{
		if (orderofControlQubit == localNum - 1)
		{
			state->Permutation(targetQubit, state->localQubits[localNum - 2]);
			orderofTargetQubit = localNum - 2;
		}
		else
		{
			state->Permutation(targetQubit, state->localQubits[localNum - 1]);
			orderofTargetQubit = localNum - 1;
		}
	}
	else if (!controlQubitLocal && targetQubitLocal)
	{
		state->whetherQubitNonlocal(state->qubitNum - localNum, controlQubit, orderofNonlocalControlQubit);
		if (!(((myid >> orderofNonlocalControlQubit) % 2) ^ (controlQubitType)))
			singleQubitGateOperates(localNum, orderofTargetQubit, state);
		return;
	}
	controlledGateOperates(localNum, controlQubitType, orderofControlQubit, orderofTargetQubit, state);
	return;
}


void Gate::operator()(int localNum, State* state, vector<int>& qubitSequence)
{
	vector<int> orderofQubitSequence(this->dim);
	int i = 0, j = 0, k = 0, indexInitial = localNum - 1;

	for (i = 0; i < this->dim; i++)
	{
		if (!state->whetherQubitLocal(localNum, qubitSequence[i], orderofQubitSequence[i]))
		{
			for (j = indexInitial; j >= 0; j--)
			{
				if (!searchElement(state->localQubits[j], this->dim, qubitSequence))
				{
					state->Permutation(qubitSequence[i], state->localQubits[j]);
					orderofQubitSequence[i] = j;
					indexInitial = j;
					break;
				}
			}
		}
	}
	generalQubitGateOperates(localNum, orderofQubitSequence, state);
	return;
}