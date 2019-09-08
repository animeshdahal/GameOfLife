
//
//
//
//
//  Created by Animesh Dahal
//
#include "mpi.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#define MASTER 0
using namespace std;


void allocate_space(int* &array,int cols){ //allocates memory for dynamically initialized arrays
    array = new int[cols];
}

int main(int argc, char* argv[])
{
    int comm_size, num_col, rank;
    
    MPI_Status Stat;
    MPI_Request send_request, recv_request;
    
    MPI_Init(&argc,&argv); //initializing MPI
    
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    num_col = atoi(argv[1]);    //input for square matrix
    int generations = atoi(argv[2]); //input
    
    int workload = num_col/comm_size;
    
    int my_start = rank * workload;
    int my_end = (rank == (comm_size-1)) ? num_col : my_start + workload;
    int my_partition = my_end - my_start;
    
    int* recvSlice;
    double startTime = MPI_Wtime();
    
    if(rank == MASTER){
        recvSlice = new int[num_col*num_col];
        
        /*Initial seed to the matrixwith int values*/
        for (unsigned int i=0; i<num_col; i++){
            for(unsigned int j=0; j<num_col; j++){
                if (i % 2 == 0){
                    if ( j % 2 == 0 )
                        recvSlice[i*num_col + j] = 1;
                    else
                        recvSlice[i*num_col + j] = 0;
                }
                else{
                    if ( j % 2 == 0 )
                        recvSlice[i*num_col + j] = 0;
                    else
                        recvSlice[i*num_col + j] = 1;
                }
            }
        }
        
        for (int worker=1; worker<comm_size; worker++)
        {
            int worker_start = worker * workload;
            int worker_end = (worker == (comm_size-1)) ? num_col : worker_start + workload;
            int worker_partition = worker_end-worker_start;
        
            MPI_Send(&(recvSlice[worker_start*num_col]), worker_partition*num_col, MPI_INT, worker, 1, MPI_COMM_WORLD);    //send initial data by master
        }
        
    }
    
    if(rank!=0){
        recvSlice = new int[my_partition * num_col];
        MPI_Recv(recvSlice, my_partition*num_col, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &Stat); //receive initial data by all process
    }
    
    /*dynamic allocation of array to pass information beteween processes*/
    int* todown;
    int* toup;
    int* tofirst;
    int* tolast;
    int* fromdown;
    int* fromup;
    int* fromlast;
    int* fromfirst;
    allocate_space(todown,num_col);
    allocate_space(toup,num_col);
    allocate_space(tofirst,num_col);
    allocate_space(tolast,num_col);
    allocate_space(fromdown,num_col);
    allocate_space(fromup,num_col);
    allocate_space(fromlast,num_col);
    allocate_space(fromfirst,num_col);
    
    
    for (int g=0; g<generations; g++) //generations loop begins
    {
        
        int prevProcess = (my_start > 0) ? rank - 1 : -1;
        int nextProcess = (my_end < (int) num_col) ? rank + 1 : -1;
        if (nextProcess != -1) {  //test if it is last row
            if (prevProcess != -1) { //test if it is first row
                
                MPI_Sendrecv(&(todown[0]),num_col,MPI_INT,nextProcess,1,&(fromup[0]),num_col,MPI_INT,prevProcess,1,MPI_COMM_WORLD, &Stat); //sending and receiving information from adjacent process
                
            }else{ //first row sends to next process, receives from last process
               
                MPI_Send(&(todown[0]), num_col, MPI_INT, nextProcess, 1, MPI_COMM_WORLD); //
                
                MPI_Recv(&(fromlast[0]), num_col, MPI_INT, comm_size-1, 1, MPI_COMM_WORLD, &Stat);
            }
        }else{
            if (prevProcess != -1) { //last row receives from previous process, sends to first process
                
                MPI_Recv(&(fromup[0]), num_col, MPI_INT, prevProcess, 1, MPI_COMM_WORLD, &Stat);

                MPI_Send(&(tofirst[0]), num_col, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
              
            }
        }

        for (int k=0; k<my_partition; k++) //for each row
        {
            for (int i=0; i<num_col; i++) //for each column
            {  //implementation of game of life logic
                if(k == 0 && i == 0){
                    int sum = recvSlice[k*num_col+i+1] + fromdown[i] + fromdown[i+1] + recvSlice[k*num_col+num_col-1] + fromlast[i] + fromlast[i+1] + fromlast[num_col-1] + fromdown[num_col-1] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k == 0 && i == num_col - 1 ){
                    int sum = recvSlice[0] + recvSlice[k*num_col+i-1] + fromdown[i-1] + fromdown[i] + fromlast[i] + fromlast[0] + fromlast[i-1] + fromdown[0] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k == num_col-1 && i == 0 ){
                    int sum = recvSlice[k*num_col+i+1] + fromup[i] + fromup[i+1] + fromfirst[i] + recvSlice[k*num_col+num_col-1] + fromfirst[num_col-1] + fromup[num_col-1] + fromfirst[i+1] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k == num_col-1 && i == num_col-1 ){
                    int sum = fromup[i] + fromup[i-1] + recvSlice[k*num_col+i-1] + fromfirst[i] + recvSlice[k*num_col+0] + fromfirst[0] + fromup[0] + fromfirst[i-1] ;
                    if ( sum < 2 || sum > 3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k == 0 && i != 0 && i != num_col-1){
                    int sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + fromdown[i] + fromdown[i-1] + fromdown[i+1] + fromlast[i] + fromlast[i-1] + fromlast[i+1] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k == num_col-1 && i != 0 && i != num_col-1){
                    int sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + fromup[i] + fromup[i-1] + fromup[i+1] + fromfirst[i] + fromfirst[i-1] + fromfirst[i+1] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(i == 0 && k != 0 && k != num_col-1){
                    int sum = fromup[i] + fromdown[i] + recvSlice[k*num_col+i+1] + fromup[i+1] + fromdown[i+1] + recvSlice[k*num_col+num_col-1] + fromup[num_col-1] + fromdown[num_col-1] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(i == num_col-1 && k != 0 && k != num_col-1){
                    int sum = fromup[i] + fromdown[i] + recvSlice[k*num_col+i-1] + fromup[i-1] + fromdown[i-1] + recvSlice[k*num_col+0] + fromup[0] + fromdown[0] ;
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
                else if(k!=0 && i!=0 && k!=my_partition-1 && i!=num_col-1){
                    int sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + recvSlice[(k-1)*num_col+i] + recvSlice[(k+1)*num_col+i] + recvSlice[(k-1)*num_col+i-1] + recvSlice[(k+1)*num_col+i-1] + recvSlice[(k-1)*num_col+i+1] + recvSlice[(k+1)*num_col+i+1];
                    if ( sum < 2 || sum >3 )
                        recvSlice[k*num_col+i] = 0;
                    else if(sum == 3)
                        recvSlice[k*num_col+i] = 1;
                }
            }
        }
        
    }
    
    
    if (rank==MASTER) //master gathers the data
    {
        for (int worker=1; worker<comm_size; worker++)
        {
            int worker_start = worker * workload;
            int worker_end = (worker == (comm_size-1)) ? num_col : worker_start + workload;
            int worker_partition = worker_end-worker_start;
        
            MPI_Recv(&(recvSlice[worker_start*num_col]), worker_partition*num_col, MPI_INT, worker, 1, MPI_COMM_WORLD, &Stat); //receive all
            
        }
    }
    else{ //process except master send respective data to master
        MPI_Send(recvSlice, num_col*my_partition, MPI_INT, MASTER,1, MPI_COMM_WORLD);
    }
    
    
    delete[] recvSlice; //delocate memory
    
    delete[] todown;
    delete[] toup;
    delete[] tofirst;
    delete[] tolast;
    delete[] fromdown;
    delete[] fromup;
    delete[] fromlast;
    delete[] fromfirst;
    
    double endTime = MPI_Wtime();
    double timeElapsed = endTime - startTime;
    if (rank == 0) { //master prints the elapsed time
        cout<<"Time Elapsed = "<<timeElapsed<<endl<<"============="<<endl<<endl;
        
    }
    
    MPI_Finalize();
}
