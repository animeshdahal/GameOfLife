//
//
//
//
//  Created by Animesh Dahal
//


#include <iostream>
#include <omp.h>

using namespace std;

int main(int argc, const char * argv[]) {
    int num_col = atoi(argv[1]); //taking matrix size from user input
    int thread_count = atoi(argv[3]); //taking number of cores from user input
    int generations = atoi(argv[2]); //taking generations input
    double** live_dead = new double*[num_col];
    for (unsigned int i=0; i<num_col; i++){
        live_dead[i] = new double[num_col];
    }
    /*Initial seed to the matrixwith double values*/
    for (unsigned int i=0; i<num_col; i++){
        for(unsigned int j=0; j<num_col; j++){
            if (i % 2 == 0){
                if ( j % 2 == 0 )
                    live_dead[i][j] = 1;
                else
                    live_dead[i][j] = 0;
            }
            else{
                if ( j % 2 == 0 )
                    live_dead[i][j] = 0;
                else
                    live_dead[i][j] = 1;
            }
        }
    }
    int partition =  num_col / thread_count; //creating partition size to partition the matrix row-wise for parallel multi-processing
    int count = 0;
    
    double parallel_startTime = omp_get_wtime(); //gets the start time
    while ( count < generations ){
#pragma omp parallel num_threads(thread_count) //parallel section starts here
    {
        int my_rank = omp_get_thread_num();
        int my_start = my_rank * partition;
        int my_end =  (my_rank == (thread_count - 1)) ? num_col : my_start + partition; //implements row-wise partition of the matrices
            for ( int k=my_start; k<my_end ; k++){
                for( int i=0; i<num_col; i++ ){
                    if(k == 0 && i == 0){
                        double sum = live_dead[k][i+1] + live_dead[k+1][i] + live_dead[k+1][i+1] + live_dead[k][num_col-1] + live_dead[num_col-1][i] + live_dead[num_col-1][i+1] + live_dead[num_col-1][num_col-1] + live_dead[k+1][num_col-1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(k == 0 && i == num_col - 1 ){
                        double sum = live_dead[k][0] + live_dead[k][i-1] + live_dead[k+1][i-1] + live_dead[k+1][i] + live_dead[num_col-1][i] + live_dead[num_col-1][0] + live_dead[num_col-1][i-1] + live_dead[k+1][0] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(k == num_col-1 && i == 0 ){
                        double sum = live_dead[k][i+1] + live_dead[k-1][i] + live_dead[k-1][i+1] + live_dead[0][i] + live_dead[k][num_col-1] + live_dead[0][num_col-1] + live_dead[k-1][num_col-1] + live_dead[0][i+1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(k == num_col-1 && i == num_col-1 ){
                        double sum = live_dead[k-1][i] + live_dead[k-1][i-1] + live_dead[k][i-1] + live_dead[0][i] + live_dead[k][0] + live_dead[0][0] + live_dead[k-1][0] + live_dead[0][i-1] ;
                        if ( sum < 2 || sum > 3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(k == 0 && i != 0 && i != num_col-1){
                        double sum = live_dead[k][i-1] + live_dead[k][i+1] + live_dead[k+1][i] + live_dead[k+1][i-1] + live_dead[k+1][i+1] + live_dead[num_col-1][i] + live_dead[num_col-1][i-1] + live_dead[num_col-1][i+1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(k == num_col-1 && i != 0 && i != num_col-1){
                        double sum = live_dead[k][i-1] + live_dead[k][i+1] + live_dead[k-1][i] + live_dead[k-1][i-1] + live_dead[k-1][i+1] + live_dead[0][i] + live_dead[0][i-1] + live_dead[0][i+1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(i == 0 && k != 0 && k != num_col-1){
                        double sum = live_dead[k-1][i] + live_dead[k+1][i] + live_dead[k][i+1] + live_dead[k-1][i+1] + live_dead[k+1][i+1] + live_dead[k][num_col-1] + live_dead[k-1][num_col-1] + live_dead[k+1][num_col-1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else if(i == num_col-1 && k != 0 && k != num_col-1){
                        double sum = live_dead[k-1][i] + live_dead[k+1][i] + live_dead[k][i-1] + live_dead[k-1][i-1] + live_dead[k+1][i-1] + live_dead[k][0] + live_dead[k-1][0] + live_dead[k+1][0] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                    else{
                        double sum = live_dead[k][i-1] + live_dead[k][i+1] + live_dead[k-1][i] + live_dead[k+1][i] + live_dead[k-1][i-1] + live_dead[k+1][i-1] + live_dead[k-1][i+1] + live_dead[k+1][i+1] ;
                        if ( sum < 2 || sum >3 )
                            live_dead[k][i] = 0;
                        else if(sum == 3)
                            live_dead[k][i] = 1;
                    }
                }
            }
        }
        count++;
        }
    
    double parallel_endTime = omp_get_wtime(); //gets the end time
    double parallel_timeElapsed = parallel_endTime - parallel_startTime; //computes total time taken for the program to perform matrix multiplication
    
    cout<<"Time Elapsed = "<<parallel_timeElapsed<<endl; //outputs the result*/
    
    for (int i=0; i<num_col; i++){
        delete[] live_dead[i];
    }
    
    delete[] live_dead;
    return 0;
    
}
