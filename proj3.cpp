/**
 * Prirazeni poradi preorder vrcholum (PRL project 3)
 * Filip Stastny (xstast24)
 */
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <algorithm>

#define BUFSIZE 128
#define TAG 0
#define NUMBERS_FILENAME "numbers"

using namespace std;

// Helper method to print vector
void printVector(vector<int> v) {
  cout << v.at(0);
  for (int i=1; i<v.size(); i++) {
      cout << " " << v.at(i);
    }
    cout << endl;
}

// Helper method to print array
void printArray(int arr[], int size) {
  cout << arr[0];
  for (int i=1; i < size; i++) {
      cout << " " << arr[i];
    }
    cout << endl;
}

// get index of given value from given array
int getElementIndex(int array[], int size, int value)
{ 
  int index = 0;
  while ( index < size && array[index] != value ) ++index;
  return ( index == size ? -1 : index );
}

int main(int argc, char *argv[])
{
  char idstr[32];
  char buff[BUFSIZE];
  MPI_Status stat;  
  int i;
  int procCount;
  int myNumberCount;
  int myid;
  string input;
  vector<int> adjacents;
  int from; // starting node
  int to; // target node
  int inverse; // inverse edge to this one
  int nextEdgeEtour; // the next Etour edge
  int weight; // weight of edge in pre-order assignment algorithm, 1=forward edge, 0=backward edge
  int suffixsum;

  MPI_Init(&argc,&argv); // initialize MPI
  MPI_Comm_size(MPI_COMM_WORLD,&procCount); // get number of running processors
  MPI_Comm_rank(MPI_COMM_WORLD,&myid); // get ID of this processor
  input = argv[argc-1]; // get the input tree string
  
  // only root given - return root
  if(procCount == 1) {
    cout << input;
    MPI_Finalize(); 
    return 0;
  }

  // get neighours for this edge
  if (myid % 2 == 0) {
    from = myid/4;
    to = myid/2+1;
  } else {
    from = myid/2+1;
    to = myid/4;
  }
  
  // create adjacency list
  int toParent = (to == 0 ? -1 : ((to-1)*2 + 1)); // root node does not have a parent
  int fromParent = (to == 0 ? -1 : ((to-1)*2)); // root node does not have a parent
  int toLeftChild = ((to+1)*4-2 > procCount) ? -1 : to*4; // test if the 
  int fromLeftChild = ((to+1)*4-2 > procCount) ? -1 : to*4 + 1;
  int toRightChild = ((to+1)*4 > procCount) ? -1 : to*4 + 2;
  int fromRightChild = ((to+1)*4 > procCount) ? -1 : to*4 + 3;
  if (toParent != -1) {adjacents.push_back(toParent);}
  if (fromParent != -1) {adjacents.push_back(fromParent);}
  if (toLeftChild != -1) {adjacents.push_back(toLeftChild);}
  if (fromLeftChild != -1) {adjacents.push_back(fromLeftChild);}
  if (toRightChild != -1) {adjacents.push_back(toRightChild);}
  if (fromRightChild != -1) {adjacents.push_back(fromRightChild);}

  // get inverse edge for this edge
  if (myid % 2 == 0) {
    inverse = myid + 1;
  } else {
    inverse = myid - 1;
  }

  // get next etour edge for this edge
  int adjacentsLen = adjacents.size();
  for (i=0; i < adjacentsLen; i = i + 2) { // iterate by i+2, cos each node has 2 edges
    if (adjacents[i] == inverse) {
      if (i >= (adjacentsLen-2)) { // next == null
        // "first item of adj. list of vertex v" (podle prednasky vezmu prvni, ale nemel by to byt jen predchozi? todo: ujasnit si)
        nextEdgeEtour = adjacents[0];
      } else { // next existuje
        nextEdgeEtour = adjacents[i+2];
      }
      break;
    }
  }

  // create Euler's Path array
  int etour[procCount];
  etour[myid] = nextEdgeEtour;
  // all processors building the etour array together
  // MPI_Allgather man page "The jth block of data sent from each process is received by every process and placed in the jth block of the buffer recvbuf"
  MPI_Allgather(&etour[myid], 1, MPI_INT, etour, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // wait for all processors to get here
  
  // get weight of this edge
  if (myid % 2 == 0) {
    weight = 1;
  } else {
    weight = 0;
  }

  // count suffix sums - linear alorithm, not parallel
  int lastEdgeIndex = getElementIndex(etour, procCount, getElementIndex(etour, procCount, 0));
  int suffixFromSender;
  if(getElementIndex(etour, procCount, myid) == etour[lastEdgeIndex]){ // root element
    MPI_Recv(&suffixFromSender, 1, MPI_INT, etour[myid], 0, MPI_COMM_WORLD, &stat);
    suffixsum = suffixFromSender + weight;
  } else if(myid == etour[lastEdgeIndex]){ // last element of etour
    suffixsum = weight;
    MPI_Send(&suffixsum, 1, MPI_INT, getElementIndex(etour, procCount, myid), 0, MPI_COMM_WORLD);
  } else {
    MPI_Recv(&suffixFromSender, 1, MPI_INT, etour[myid], 0, MPI_COMM_WORLD, &stat);
    suffixsum = suffixFromSender + weight;
   MPI_Send(&suffixsum, 1, MPI_INT, getElementIndex(etour, procCount, myid), 0, MPI_COMM_WORLD);
  }
  suffixsum = suffixsum * weight;
  MPI_Barrier(MPI_COMM_WORLD); // wait for all processors to get here

  // gather suffixes and print the tree
  if (myid == 0) {
    int tmpSuffix = 0;
    int tmpNode = 0;
    int results[procCount];
    int nodes[procCount];
    results[0] = suffixsum;
    nodes[0] = to;
    for(i=1; i < procCount; i++) {
      MPI_Recv(&tmpSuffix, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &stat);
      MPI_Recv(&tmpNode, 1, MPI_INT, i, TAG, MPI_COMM_WORLD, &stat);
      results[i] = tmpSuffix;
      nodes[i] = tmpNode;
    }

    // print nodes by value (first the greatest ones)
    cout << input[0];
    for (i=0; i<input.length()-1; i++) {
      // get max element
      int tmpMax = -1;
      int tmpIndex = -1;
      for (int j=0; j<procCount; j++) {
        if (results[j] > tmpMax) {
          tmpMax = results[j];
          tmpIndex = j;
        }
      }
      results[tmpIndex] = -1; // delete the taken element from results
      int resultIndex = nodes[tmpIndex];
      cout << input[resultIndex];
    }
  } else {
    MPI_Send(&suffixsum, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
    MPI_Send(&to, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
  }

  MPI_Finalize(); 
  return 0;
}