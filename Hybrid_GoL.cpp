
//
//
//
//
//  Created by Animesh Dahal
//
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <fstream>
#include <sstream>
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

    string infilename;
    string outfilename;
    infilename = argv[1];
    int generations = atoi(argv[2]);
    int thread_count = atoi(argv[3]);
    outfilename = argv[4];
    
    int* recvSlice;
    int* tempSlice;
    string line;
    int x;
    int y;
    int rowIP;
    int colIP;
    int count = 0;
    
    
    ifstream inFile;
    ofstream outFile;
    
    if(rank == MASTER){//Master opens the input file
        int row = 0;
        int col = 0;
        inFile.open(infilename);
        
        if(inFile.fail()){
            cerr << "File cannot be opened"<<endl;
            exit(1);
        }
        
        while(inFile.good()){ //Master streams the data from the input file
            if(count == 0){
                getline(inFile, line);
                istringstream streamA(line);
                streamA >> x >> y;
                rowIP = x;
                colIP = y;
                recvSlice = new int[rowIP*colIP];
                tempSlice = new int[rowIP*colIP];
            }
            getline(inFile, line);
            istringstream streamA(line);
            col = 0;
            while(streamA >> x){
                recvSlice[row * colIP + col] = x;
                col++;
            }
            row++;
            
            count++;
        }
        for (unsigned int i=0; i<rowIP*colIP; i++)
            tempSlice[i] = recvSlice[i];
        
        for (int worker=1; worker<comm_size; worker++){ //master sends row information to all processes
            MPI_Send(&colIP, 4, MPI_INT, worker, 0, MPI_COMM_WORLD);
        }
    }
    
    else{//all processes receive row information
        int num_coloumn;
        MPI_Recv(&num_coloumn, 4, MPI_INT, MASTER, 0, MPI_COMM_WORLD, &Stat);
        num_col = num_coloumn;
    }
    if (rank == MASTER)
        num_col = colIP;
    
    if(rank == MASTER){//opens and writes coloumn and row information to output file
        outFile.open(outfilename);

        outFile << rowIP <<" "<<colIP<<endl;
    }
    
    double startTime = MPI_Wtime();
    
    int workload = num_col/comm_size;
    int my_start = rank * workload;
    int my_end = (rank == (comm_size-1)) ? num_col : my_start + workload;
    int my_partition = my_end - my_start;
    
    if(rank == MASTER){
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
        tempSlice = new int[my_partition * num_col];
        MPI_Recv(recvSlice, my_partition*num_col, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &Stat); //receive initial data by all processes
        for(int i = 0; i<my_partition*num_col; i++){
            tempSlice[i] = recvSlice[i];
        }
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
                for(int m = 0; m<num_col; m++){
                    toup[m] = recvSlice[m];
                    todown[m] = recvSlice[ (my_partition-1) * num_col + m ];
                }
                
                MPI_Sendrecv(&(todown[0]),num_col,MPI_INT,nextProcess,1,&(fromup[0]),num_col,MPI_INT,prevProcess,1,MPI_COMM_WORLD, &Stat); //sending and receiving information from adjacent process
                
                MPI_Recv(&(fromdown[0]), num_col, MPI_INT, rank+1, 1, MPI_COMM_WORLD, &Stat);
                
                MPI_Send(&(toup[0]),num_col,MPI_INT,rank-1,1,MPI_COMM_WORLD);
                
            }else{ //first row sends to and receive from next and last processes
                for(int m = 0; m<num_col; m++){
                    tolast[m] = recvSlice[m];
                    todown[m] = recvSlice[ (my_partition-1) * num_col + m ];
                }
                MPI_Send(&(todown[0]), num_col, MPI_INT, nextProcess, 1, MPI_COMM_WORLD); //
                
                MPI_Send(&(tolast[0]), num_col, MPI_INT, comm_size-1, 1, MPI_COMM_WORLD);
                
                MPI_Recv(&(fromlast[0]), num_col, MPI_INT, comm_size-1, 1, MPI_COMM_WORLD, &Stat);
                
                MPI_Recv(&(fromdown[0]), num_col, MPI_INT, rank + 1, 1, MPI_COMM_WORLD, &Stat);
            }
        }else{
            if (prevProcess != -1) { //last row sends to and receives from previous and firt process
                for(int m = 0; m<num_col; m++){
                    tofirst[m] = recvSlice[ (my_partition-1) * num_col + m ];
                    toup[m] = recvSlice[m];
                }
                
                MPI_Recv(&(fromup[0]), num_col, MPI_INT, prevProcess, 1, MPI_COMM_WORLD, &Stat);
                
                MPI_Recv(&(fromfirst[0]), num_col, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &Stat);
                
                MPI_Send(&(tofirst[0]), num_col, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
                
                MPI_Send(&(toup[0]), num_col, MPI_INT, rank-1 , 1, MPI_COMM_WORLD);
            }
        }
#pragma omp parallel num_threads(thread_count)//multi-threading begins
        {
            int my_workload_omp = workload / thread_count;
            int my_rank_omp = omp_get_thread_num();
            int my_start_omp = my_rank_omp * my_workload_omp;
            
            int my_end_omp = (my_rank_omp == (thread_count-1)) ? my_partition : my_start_omp + my_workload_omp;
            
            int partition =  my_end_omp - my_start_omp ;
            
            for (int k=my_start_omp; k<my_end_omp; k++) //for each row
            {
                for (int i=0; i<num_col; i++) //for each column
                {  //implementation of game of life logic
                    if(k == 0 && i == 0){
                        int sum;
                        if(rank == MASTER){
                            sum = recvSlice[k*num_col+i+1] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col+i+1] + recvSlice[k*num_col+num_col-1] + fromlast[i] + fromlast[i+1] + fromlast[num_col-1] + recvSlice[(k+1)*num_col+num_col-1] ;
                        }
                        else{
                            sum = recvSlice[k*num_col+i+1] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col+i+1] + recvSlice[k*num_col+num_col-1] + fromup[i] + fromup[i+1] + fromup[num_col-1] + recvSlice[(k+1)*num_col+num_col-1] ;
                        }
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else if( k == 0 && i == num_col - 1 ){
                        int sum;
                        if(rank == MASTER)
                            sum = recvSlice[0] + recvSlice[k*num_col+i-1] + recvSlice[(k+1)*num_col+i-1] + recvSlice[(k+1)*num_col+i] + fromlast[i] + fromlast[0] + fromlast[i-1] + recvSlice[(k+1)*num_col] ;
                        else
                            sum = recvSlice[0] + recvSlice[k*num_col+i-1] + recvSlice[(k+1)*num_col+i-1] + recvSlice[(k+1)*num_col+i] + fromup[i] + fromup[0] + fromup[i-1] + recvSlice[(k+1)*num_col] ;
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else if( i == 0 && k == my_partition - 1){
                        int sum;
                        if (rank == comm_size-1)
                            sum = recvSlice[k*num_col+i+1] + recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i+1] + fromfirst[i] + recvSlice[k*num_col+num_col-1] + fromfirst[num_col-1] + recvSlice[(k-1)*num_col+num_col-1] + fromfirst[i+1] ;
                        else
                            sum = recvSlice[k*num_col+i+1] + recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i+1] + fromdown[i] + recvSlice[k*num_col+num_col-1] + fromdown[num_col-1] + recvSlice[(k-1)*num_col+num_col-1] + fromdown[i+1] ;
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else if( i == num_col-1 &&  k == my_partition -1 ){
                        int sum;
                        if (rank == comm_size-1)
                            sum = recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i-1] + recvSlice[k*num_col+i-1] + fromfirst[i] + recvSlice[k*num_col+0] + fromfirst[0] + recvSlice[(k-1)*num_col] + fromfirst[i-1] ;
                        else
                            sum = recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i-1] + recvSlice[k*num_col+i-1] + fromdown[i] + recvSlice[k*num_col+0] + fromdown[0] + recvSlice[(k-1)*num_col] + fromdown[i-1] ;
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else if(k == 0 && i > 0 && i < num_col-1){
                        int sum;
                        if(rank == MASTER){
                            sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col+i-1] + recvSlice[(k+1)*num_col+i+1] + fromlast[i] + fromlast[i-1] + fromlast[i+1] ;
                        }
                        else
                            sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col+i-1] + recvSlice[(k+1)*num_col+i+1] + fromup[i] + fromup[i-1] + fromup[i+1] ;
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else if( k == my_partition-1 && i > 0 && i < num_col-1){
                        int sum;
                        if (rank == comm_size -1)
                            sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i-1] + recvSlice[(k-1)*num_col+i+1] + fromfirst[i] + fromfirst[i-1] + fromfirst[i+1] ;
                        else
                            sum = recvSlice[k*num_col+i-1] + recvSlice[k*num_col+i+1] + recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i-1] + recvSlice[(k-1)*num_col+i+1] + fromdown[i] + fromdown[i-1] + fromdown[i+1] ;
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                    
                    else {
                        int sum;
                        if(i==0 && k < my_partition && k > 0){
                            sum = recvSlice[(k-1)*num_col+i+1] + recvSlice[(k-1)*num_col+num_col-1] + + recvSlice[(k-1)*num_col+0] + recvSlice[(k+1)*num_col+i+1] + recvSlice[(k+1)*num_col+0] + recvSlice[(k+1)*num_col+num_col-1] + recvSlice[k*num_col+num_col-1] + recvSlice[k*num_col+i+1]   ;
                        }
                        else if(i==num_col-1 && my_partition > k && k > 0 ){
                            sum = recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col] + recvSlice[(k-1)*num_col+i-1] + recvSlice[k*num_col+i-1] +  recvSlice[k*num_col] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col] + recvSlice[(k+1)*num_col+i-1];
                        }
                        else if(my_partition-1 > k && k > 0) {
                            sum = recvSlice[(k-1)*num_col+i] + recvSlice[(k-1)*num_col+i+1] + recvSlice[(k-1)*num_col+i-1] + recvSlice[k*num_col+i-1] +  recvSlice[k*num_col+i+1] + recvSlice[(k+1)*num_col+i] + recvSlice[(k+1)*num_col+i+1] + recvSlice[(k+1)*num_col+i-1];
                        }
                        if ( sum < 2 || sum >3 )
                            tempSlice[k*num_col+i] = 0;
                        else if(sum == 3)
                            tempSlice[k*num_col+i] = 1;
                    }
                }
            }
        }
        for(int i = 0; i<my_partition*num_col; i++){
            recvSlice[i] = tempSlice[i];
        }
    }
        
        if (rank==MASTER) //master gathers the data
        {
            for (int worker=1; worker<comm_size; worker++)
            {
                int worker_start = worker * workload;
                int worker_end = (worker == (comm_size-1)) ? num_col : worker_start + workload;
                int worker_partition = worker_end - worker_start;
                
                MPI_Recv(&(recvSlice[worker_start*num_col]), worker_partition*num_col, MPI_INT, worker, 1, MPI_COMM_WORLD, &Stat); //receive all
                
            }
        }
        else{ //process except master send respective data to master
            MPI_Send(recvSlice, num_col*my_partition, MPI_INT, MASTER,1, MPI_COMM_WORLD);
        }
    
    double endTime = MPI_Wtime();
    double timeElapsed = endTime - startTime;
    if (rank == 0) { //master prints the elapsed time
        cout<<"Time Elapsed = "<<timeElapsed<<endl<<"============="<<endl<<endl;
        
    }
    
    
    
        if (rank==MASTER){//master writes data into outputfile
            for (unsigned int i=0; i<num_col; i++){
                for(unsigned int j=0; j<num_col; j++){
                    outFile << recvSlice[i*num_col + j] << " ";
                }
                outFile << endl;
            }
            inFile.close();
            outFile.close();
        }
        
        delete[] recvSlice; //delocate memory
        delete[] tempSlice;
        
        delete[] todown;
        delete[] toup;
        delete[] tofirst;
        delete[] tolast;
        delete[] fromdown;
        delete[] fromup;
        delete[] fromlast;
        delete[] fromfirst;
        
    
        
        MPI_Finalize();
    }
