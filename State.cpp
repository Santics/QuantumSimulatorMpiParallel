#include "State.h"
#define QUBITMAXSIZE  50

//***********************
//developer based functions

//judge whether a specific qubit is local
bool State::whetherQubitLocal(int localNum, int qubit, int& orderofAimQubit)
{
	int i = 0;
	bool qbitLocal = false;

	orderofAimQubit = -1;
	for (i = 0; i < localNum; i++)
	{
		if (qubit == this->localQubits[i])
		{
			qbitLocal = true;
			orderofAimQubit = i;
			break;
		}
	}

	return qbitLocal;
}


//judge whether a specific qubit is nonlocal
bool State::whetherQubitNonlocal(int nonlocalNum, int qubit, int& orderofAimQubit)
{
	int i = 0;
	bool qbitnonLocal = false;

	orderofAimQubit = -1;
	for (i = 0; i < nonlocalNum; i++)
	{
		if (qubit == this->nonLocalQubits[i])
		{
			qbitnonLocal = true;
			orderofAimQubit = i;
			break;
		}
	}

	return qbitnonLocal;
}


void State::LocateQubit(int nonlocalQubit, int localQubit, int nonlocalQubitNum, int localQubitNum, int& nonlocalIndex, int& localIndex)
{
	int i, j;

	for (i = 0; i < nonlocalQubitNum; i++)
	{
		if (nonlocalQubit == this->nonLocalQubits[i])
			break;
	}
	nonlocalIndex = i;

	for (j = 0; j < localQubitNum; j++)
	{
		if (localQubit == this->localQubits[j])
			break;
	}
	localIndex = j;

	return;

}

int State::qubitFind(vector<int> qubits, int  targetQubit)
{
	int i = 0;

	for (i = 0; i < qubits.size(); i++)
	{
		if (qubits[i] == targetQubit)
			return i;
	}

	return -1;
}




//***********************
//user based functions
State::State(int qubitNum)
{

	int dim = 0;
	int numprocs = 0;
	int nonlocalQubitNum = 0;
	int count = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	nonlocalQubitNum = log(numprocs) / log(2);
	dim = pow(2, qubitNum - nonlocalQubitNum);


	this->qubitNum = qubitNum;
	this->nonLocalQubits = (int*)operator new(nonlocalQubitNum * sizeof(int));
	this->localQubits = (int*)operator new((qubitNum - nonlocalQubitNum) * sizeof(int));
	this->stateVector = VectorXcd::Zero(dim);

	for (count = 0; count < qubitNum - nonlocalQubitNum; count++)
		localQubits[count] = count;
	for (count = qubitNum - nonlocalQubitNum; count < qubitNum; count++)
		nonLocalQubits[count - (qubitNum - nonlocalQubitNum)] = count;


	return;
}

void State::Initialize(std::string initialQubitString)
{
	int numprocs = 0;
	int myid = 0;
	int nonlocalQubitNum = 0, localQubitNum = 0;
	int nonLocalIndex = 0, localIndex = 0;


	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);


	nonlocalQubitNum = log(numprocs) / log(2);
	localQubitNum = this->qubitNum - nonlocalQubitNum;

	bitset<QUBITMAXSIZE> initialLocalQubitSequence(initialQubitString, nonlocalQubitNum, localQubitNum);
	bitset<QUBITMAXSIZE> initialNonlocalQubitSequence(initialQubitString, 0, nonlocalQubitNum);

	localIndex = initialLocalQubitSequence.to_ulong();
	nonLocalIndex = initialNonlocalQubitSequence.to_ulong();


	if (myid == nonLocalIndex)
	{
		this->stateVector(localIndex) = 1;
	}

	return;

}


void State::Permutation(int nonlocalQubit, int localQubit)
{
	int nonlocalIndex, localIndex;
	int nonlocalPartition, nonlocalPartitionElements;
	int localPartition, localPartitionElements;
	int numprocs = 0, myid = 0;
	int nonlocalQubitNum = 0, localQubitNum = 0;
	int i, j, k;
	MPI_Status status;
	VectorXcd buffer;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	if (numprocs % 2 == 0)
	{

		nonlocalQubitNum = log(numprocs) / log(2);
		localQubitNum = this->qubitNum - nonlocalQubitNum;
		LocateQubit(nonlocalQubit, localQubit, nonlocalQubitNum, localQubitNum, nonlocalIndex, localIndex);

		nonlocalPartition = pow(2, (nonlocalQubitNum - nonlocalIndex - 1));
		nonlocalPartitionElements = pow(2, nonlocalIndex + 1);
		localPartition = pow(2, (localQubitNum - localIndex - 1));
		localPartitionElements = pow(2, localIndex + 1);

		for (i = 0; i < localPartition; i++)

		{

			if (((myid / (int)pow(2, nonlocalIndex)) % 2) == 0)
			{
				MPI_Send(&this->stateVector[(i * localPartitionElements) + (localPartitionElements / 2)], localPartitionElements / 2, MPI_C_DOUBLE_COMPLEX, myid + (nonlocalPartitionElements / 2), 0, MPI_COMM_WORLD);
			}
			else
			{
				buffer = VectorXd::Zero(localPartitionElements / 2);
				MPI_Recv(buffer.data(), localPartitionElements / 2, MPI_C_DOUBLE_COMPLEX, myid - (nonlocalPartitionElements / 2), 0, MPI_COMM_WORLD, &status);
				MPI_Send(&this->stateVector.data()[i * localPartitionElements], localPartitionElements / 2, MPI_C_DOUBLE_COMPLEX, myid - (nonlocalPartitionElements / 2), 0, MPI_COMM_WORLD);
				for (k = 0, j = i * localPartitionElements; k < localPartitionElements / 2; k++, j++)
					this->stateVector(j) = buffer(k);
			}
			if (((myid / (int)pow(2, nonlocalIndex)) % 2) == 0)
			{
				MPI_Recv(&this->stateVector.data()[(i * localPartitionElements) + (localPartitionElements / 2)], localPartitionElements / 2, MPI_C_DOUBLE_COMPLEX, myid + (nonlocalPartitionElements / 2), 0, MPI_COMM_WORLD, &status);
			}
		}

		this->localQubits[localIndex] = nonlocalQubit;
		this->nonLocalQubits[nonlocalIndex] = localQubit;

		return;
	}

}


void State::Permutation(int permutationNum, int* nonlocalQubit, int* localQubit)
{

}


void State::measure(int localNum, int measureQubit, double* probability)
{
	int i, j;
	bool local;
	int orderofAimQubit;
	int myid;
	int nonlocalNum = this->qubitNum - localNum;
	int partitionNum, elementNum;
	double localProbability[2] = { 0 };
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//qubit to be measured is local
	if (whetherQubitLocal(localNum, measureQubit, orderofAimQubit))
	{
		partitionNum = pow(2, localNum - orderofAimQubit);
		elementNum = pow(2, orderofAimQubit);

		for (i = 0; i < partitionNum; i++)
		{
			if (i % 2 == 0)
			{
				for (j = 0; j < elementNum; j++)
				{
					{
						localProbability[0] += norm(this->stateVector(i * elementNum + j));
					}
				}

			}
			else
			{
				for (j = 0; j < elementNum; j++)
				{
					{
						localProbability[1] += norm(this->stateVector(i * elementNum + j));
					}
				}
			}
		}
	}

	//qubit to be measured is nonlocal
	else if (whetherQubitNonlocal(this->qubitNum - localNum, measureQubit, orderofAimQubit))
	{

		if ((int(myid / pow(2, orderofAimQubit)) % 2) == 0)
		{
			for (i = 0; i < pow(2, localNum); i++)
				localProbability[0] += norm(this->stateVector(i));
			//cout << "myid: " << myid << "  local0: " << localProbability[0] << endl;
		}
		else
		{
			for (i = 0; i < pow(2, localNum); i++)
				localProbability[1] += norm(this->stateVector(i));
			//cout << "myid: " << myid << "  local1: " << localProbability[1] << endl;
		}
	}

	MPI_Reduce(localProbability, probability, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return;
}


//void State::measure(int localNum, vector<int>& measureQubit, std::string qubits,double& probability)
//{
//	vector<int> measureLocalQubitOrder, measureLocalQubits,measureNonlocalQubitOrder, measureNonlocalQubits, localResOrder;
//	int i = 0, j = 0, k = 0, m = 0, nonlocalTargetIndex = 0, localTargetIndex = 0, localQubitIndex = 0, nonlocalQubitIndex = 0, temp = 0;
//	int totalQubitNum = measureQubit.size(), indice = 0, baseIndice = 0;
//	int numprocs, myid;
//	bool whetherMeasure;
//	double localProbability = 0;
//
//	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
//
//
//	for (i = 0; i < totalQubitNum; i++)
//	{
//
//		if (whetherQubitLocal(localNum, measureQubit[i], temp))
//		{
//			measureLocalQubitOrder.push_back(temp);
//			if (qubits[totalQubitNum-i-1] == '1')
//				localTargetIndex += pow(2, i);
//			measureLocalQubits.push_back(measureQubit[i]);
//		}
//		else if (whetherQubitNonlocal(this->qubitNum - localNum, measureQubit[i], temp))
//		{
//			measureNonlocalQubitOrder.push_back(temp);
//			if (qubits[totalQubitNum-i-1] == '1')
//				nonlocalTargetIndex += pow(2, i);
//			measureNonlocalQubits.push_back(measureQubit[i]);
//		}
//	}
//
//
//	//nonlocal part
//	if (!measureNonlocalQubitOrder.empty())
//	{
//		//calculate critical parameter:nonlocalQubitIndex for each node
//		for (i = 0; i < measureNonlocalQubitOrder.size(); i++)
//		{
//			if(myid & (int)pow(2, measureNonlocalQubitOrder[i]))
//				nonlocalQubitIndex += pow(2, qubitFind(measureQubit, measureNonlocalQubits[i]));
//		}
//	}
//	else
//		nonlocalQubitIndex = -1;
//
//
//	//localpart
//	if (!measureLocalQubitOrder.empty() && (myid == nonlocalQubitIndex ||  (nonlocalQubitIndex == -1)))
//	{
//		//find qubits which is local but not to be measured
//		for (i = 0; i < localNum; i++)
//		{
//			whetherMeasure = false;
//			for (j = 0; j < measureLocalQubitOrder.size(); j++)
//			{
//				if (i == measureLocalQubitOrder[j])
//				{
//					whetherMeasure = true;
//					break;
//				}
//			}
//
//			if (!whetherMeasure)
//				localResOrder.push_back(i);
//		}
//
//
//		//calculate local probability
//		for (i = 0; i < pow(2, localNum - measureLocalQubitOrder.size()); i++)
//		{
//			baseIndice = 0;
//			for (j = 1, k = 0; k < localNum - measureLocalQubitOrder.size(); j = j << 1, k++)
//			{
//				if (i & j)
//					baseIndice += pow(2, localResOrder[k]);
//			}
//
//			indice = baseIndice;
//
//			for (j = 0; j < measureLocalQubitOrder.size(); j++)
//			{
//				indice = baseIndice;
//				//if ( localTargetIndex & (int)pow(2,(find(measureQubit.begin(),measureQubit.end(),localQubits[j])-measureQubit.begin())))
//				if(localTargetIndex & (int)pow(2,qubitFind(measureQubit,measureLocalQubits[j])))
//				{
//					indice += pow(2, measureLocalQubitOrder[j]);
//				}
//			}
//
//			localProbability += norm(this->stateVector(indice));
//		}
//	}
//
//	else if(measureLocalQubitOrder.empty() && myid == nonlocalQubitIndex)
//	{
//		for (i = 0; i < localNum; i++)
//		{
//			localProbability += norm(this->stateVector[i]);
//		}
//	}
//
//
//	MPI_Reduce(&localProbability, &probability, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//
//	return;
//}


double* State::measure(int localNum, vector<int>& measureQubit)
{
	vector<int> measureLocalQubitOrder, measureLocalQubits, measureNonlocalQubitOrder, measureNonlocalQubits, localResOrder;
	int i = 0, j = 0, k = 0, m = 0, temp = 0, localQubitIndex = 0, nonlocalQubitIndex = 0, baseIndice = 0, indice = 0;
	int numprocs, myid;
	double* localProbability = NULL, * probability = NULL;
	bool whetherMeasure;
	int totalQubitNum = measureQubit.size();

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	localProbability = new double[(int)pow(2, totalQubitNum)]();

	//divide local & nonlocal qubits
	//store their qubitorder
	for (i = 0; i < measureQubit.size(); i++)
	{

		if (whetherQubitLocal(localNum, measureQubit[i], temp))
		{
			measureLocalQubitOrder.push_back(temp);
			measureLocalQubits.push_back(measureQubit[i]);
		}
		else if (whetherQubitNonlocal(this->qubitNum - localNum, measureQubit[i], temp))
		{
			measureNonlocalQubitOrder.push_back(temp);
			measureNonlocalQubits.push_back(measureQubit[i]);
		}
	}

	//nonlocal part
	if (!measureNonlocalQubitOrder.empty())
	{
		//calculate critical parameter:nonlocalQubitIndex for each node
		for (i = 0; i < measureNonlocalQubitOrder.size(); i++)
		{
			if (myid & ((int)pow(2, measureNonlocalQubitOrder[i])))
			{
				nonlocalQubitIndex += pow(2, qubitFind(measureQubit, measureNonlocalQubits[i]));
			}
		}
	}
	else
		nonlocalQubitIndex = 0;


	//localpart
	if (!measureLocalQubitOrder.empty())
	{
		//find qubits which is local but not to be measured
		for (i = 0; i < localNum; i++)
		{
			whetherMeasure = false;
			for (j = 0; j < measureLocalQubitOrder.size(); j++)
			{
				if (i == measureLocalQubitOrder[j])
				{
					whetherMeasure = true;
					break;
				}
			}
			if (!whetherMeasure)
				localResOrder.push_back(i);
		}

		//calculate local probability
		for (i = 0; i < pow(2, localNum - measureLocalQubitOrder.size()); i++)
		{
			baseIndice = 0;
			for (j = 1, k = 0; k < localNum - measureLocalQubitOrder.size(); j = j << 1, k++)
			{
				if (i & j)
					baseIndice += pow(2, localResOrder[k]);
			}

			for (j = 0; j < pow(2, measureLocalQubitOrder.size()); j++)
			{
				localQubitIndex = 0;
				indice = baseIndice;
				for (k = 1, m = 0; m < measureLocalQubitOrder.size(); k = k << 1, m++)
				{
					if (j & k)
					{
						indice += pow(2, measureLocalQubitOrder[m]);
						localQubitIndex += pow(2, qubitFind(measureQubit, measureLocalQubits[m]));
					}
				}
				localProbability[nonlocalQubitIndex + localQubitIndex] += norm(this->stateVector[indice]);
			}
		}
	}

	else
	{
		for (i = 0; i < pow(2, localNum); i++)
		{
			localProbability[nonlocalQubitIndex] += norm(this->stateVector[i]);
		}

	}

	//MPI communicate & reduce part
	if (myid == 0)
	{
		probability = new double[(int)pow(2, totalQubitNum)]();
	}

	MPI_Reduce(localProbability, probability, (int)pow(2, totalQubitNum), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return probability;
}

void State::measure(int localNum, vector<int>& measureQubit, std::string qubits, double& probability)
{
	vector<int> measureLocalQubitOrder, measureLocalQubits, measureNonlocalQubitOrder, measureNonlocalQubits, localResOrder;
	int i = 0, j = 0, k = 0, m = 0, nonlocalTargetIndex = 0, localTargetIndex = 0, localQubitIndex = 0, nonlocalQubitIndex = 0, temp = 0;
	int totalQubitNum = measureQubit.size(), indice = 0, baseIndice = 0;
	int numprocs, myid;
	bool whetherMeasure;
	double localProbability = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);


	for (i = 0; i < totalQubitNum; i++)
	{

		if (whetherQubitLocal(localNum, measureQubit[i], temp))
		{
			measureLocalQubitOrder.push_back(temp);
			if (qubits[totalQubitNum - i - 1] == '1')
				localTargetIndex += pow(2, i);
			measureLocalQubits.push_back(measureQubit[i]);
		}
		else if (whetherQubitNonlocal(this->qubitNum - localNum, measureQubit[i], temp))
		{
			measureNonlocalQubitOrder.push_back(temp);
			if (qubits[totalQubitNum - i - 1] == '1')
				nonlocalTargetIndex += pow(2, i);
			measureNonlocalQubits.push_back(measureQubit[i]);
		}
	}


	//nonlocal part
	if (!measureNonlocalQubitOrder.empty())
	{
		//calculate critical parameter:nonlocalQubitIndex for each node
		for (i = 0; i < measureNonlocalQubitOrder.size(); i++)
		{
			if (myid & (int)pow(2, measureNonlocalQubitOrder[i]))
				nonlocalQubitIndex += pow(2, qubitFind(measureQubit, measureNonlocalQubits[i]));
		}
	}
	else
		nonlocalQubitIndex = -1;


	//localpart
	if (!measureLocalQubitOrder.empty() && ((nonlocalTargetIndex == nonlocalQubitIndex) || (nonlocalQubitIndex == -1)))
	{
		//find qubits which is local but not to be measured
		for (i = 0; i < localNum; i++)
		{
			whetherMeasure = false;
			for (j = 0; j < measureLocalQubitOrder.size(); j++)
			{
				if (i == measureLocalQubitOrder[j])
				{
					whetherMeasure = true;
					break;
				}
			}

			if (!whetherMeasure)
				localResOrder.push_back(i);
		}


		//calculate local probability
		indice = 0;
		for (j = 0; j < measureLocalQubitOrder.size(); j++)
		{
			if (localTargetIndex & (int)pow(2, qubitFind(measureQubit, measureLocalQubits[j])))
			{
				indice += pow(2, measureLocalQubitOrder[j]);
			}
		}
		for (i = 0; i < pow(2, localNum - measureLocalQubitOrder.size()); i++)
		{
			baseIndice = 0;
			for (j = 1, k = 0; k < localNum - measureLocalQubitOrder.size(); j = j << 1, k++)
			{
				if (i & j)
					baseIndice += pow(2, localResOrder[k]);
			}
			localProbability += norm(this->stateVector(baseIndice + indice));
		}
	}

	else if (measureLocalQubitOrder.empty() && (nonlocalTargetIndex == nonlocalQubitIndex))
	{
		for (i = 0; i < localNum; i++)
		{
			localProbability += norm(this->stateVector[i]);
		}
	}


	MPI_Reduce(&localProbability, &probability, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	return;
}