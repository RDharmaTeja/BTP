#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>

#define dataType long double
#define mpiDType MPI::LONG_DOUBLE
#define mpiInt MPI::INT
#define mpiSend MPI::COMM_WORLD.Send
#define mpiRecv MPI::COMM_WORLD.Recv
#define mpiISend MPI::COMM_WORLD.Irecv
#define mpiIRecv MPI::COMM_WORLD.Isend

using namespace std;


/******************** supports start************************/

struct  flow_variables
{
    /* data */
    dataType press;
    dataType dens;
    dataType temp;
    dataType u;
    dataType v;
    dataType w;
    dataType mach;
};

struct grid_index
{
// structure to store grid indexes
    int x;
    int y;
};

void swap (int *a, int *b)
{
// swapping two elements
    int temp = *a;
    *a = *b;
    *b = temp;
}
/******************** supports end************************/

/*********************main class start*************************/
class domain
{
//Main class for the process, it sets all the required parameters to work with
// domain parameters, grid points, num of processes in horiz and vertical
public:
    int width; // horizontal grid points
    int height; // vertical grid points
    int horz_proc;// horizontal processes (cols)
    int vert_proc;// vertical processes (rows)
//so total processes will be (h*v + 2)
    int left_rank,top_rank,right_rank,bottom_rank,current_rank;//ranks of surrounding processes computed by Set_ranks
    grid_index grid_bottom,grid_top; //stores the index of the grid points that the process is operating on, computed by set_indexes
// grid data X
    std::vector <dataType> X;
    std::vector <dataType> Y;
// flow prop
    flow_variables init_flow;
// prototype declarations
    void Set_params(int a, int b, int c, int d);// set the parameters
    void Set_ranks(int rank,int size); // sets surrounding processes
    void Set_indexes(); // sets grid indexres
    void Send_message(dataType *A,int size,int dest, int tag);
    void Recv_message(dataType *A,int max_size, int source, int tag);
    void Send_to_surr(int n);
    void Send_to_surr2();
    void Load_grid();
    /*
    order of running
    1. Set_params()
    2. Load_grid()
    3. Set_ranks()
    4. Set_indexes()
    5.
    */
};


void domain::Load_grid()
{
// to load grid into an array
    ifstream grid;
    grid.open("bumpgrid.txt");
    int imax,jmax;
    dataType a,b;
    grid >> imax >> jmax;
    width = imax;
    height = jmax;
// loading grid points into memeber variable
    while(grid >> a >> b)
    {
        X.push_back(a);
        Y.push_back(b);
    }
    cout << X[1] << endl;
}


void domain::Send_to_surr(int n)
{
    if(current_rank > 0)
    {
//int n = current_rank*current_rank;
        int tag = 99,rec;
        if(left_rank > 0)
        {
            mpiSend(&n,1,mpiInt,left_rank,tag);
            mpiRecv(&rec,1,mpiInt,left_rank,tag);
            cout << "This is "<< current_rank << " received "<< rec << " from left "<< left_rank  << endl;
        }
        if(top_rank > 0)
        {
            mpiSend(&n,1,mpiInt,top_rank,tag);
            mpiRecv(&rec,1,mpiInt,top_rank,tag);
            cout << "This is "<< current_rank << " received "<< rec << " from top "<< top_rank  << endl;
        }
        if(right_rank > 0)
        {
            mpiRecv(&rec,1,mpiInt,right_rank,tag);
            mpiSend(&n,1,mpiInt,right_rank,tag);
            cout << "This is "<< current_rank << " received "<< rec << " from right "<< right_rank  << endl;
        }
        if(bottom_rank > 0)
        {
            mpiRecv(&rec,1,mpiInt,bottom_rank,tag);
            mpiSend(&n,1,mpiInt,bottom_rank,tag);
            cout << "This is "<< current_rank << " received "<< rec << " from bottom "<< bottom_rank  << endl;
        }

    }
}


void domain::Send_to_surr2()
{
    MPI::Request r1[8];
    int n = current_rank*current_rank;
    int buf[]= {0,0,0,0,0,0,0,0};
    int rec[4];
    int tag =1;
    if(current_rank>1)
    {

        if(current_rank%2 == 1 || 1)
        {
//cout << " current rank is"<< current_rank << endl;

            if(left_rank!= 1)
            {
                buf[0]++;
                r1[0] = mpiISend(&n,1,mpiInt,left_rank,tag);
                //MPI::COMM_WORLD.Recv(&rec,1,mpiInt,left_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< left_rank <<endl ;
            }
            if(top_rank!= 1)
            {
                buf[1]++;
                r1[1] = mpiISend(&n,1,mpiInt,top_rank,tag);
                //MPI::COMM_WORLD.Recv(&rec,1,mpiInt,top_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< top_rank  << endl;
            }
            if(bottom_rank!= 1)
            {
                buf[2]++;
                r1[2] = mpiISend(&n,1,mpiInt,bottom_rank,tag);
                //MPI::COMM_WORLD.Recv(&rec,1,mpiInt,bottom_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< bottom_rank  << endl;
            }
            if(right_rank!= 1)
            {
                buf[3]++;
                r1[3] = mpiISend(&n,1,mpiInt,right_rank,tag);
                //MPI::COMM_WORLD.Recv(&rec,1,mpiInt,right_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< right_rank <<endl;
            }

            if(left_rank!= 1)
            {
                //MPI::COMM_WORLD.Send(&n,1,mpiInt,left_rank,tag);
                buf[4]++;
                r1[4] = mpiIRecv(rec,1,mpiInt,left_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< left_rank <<endl ;
            }
            if(top_rank!= 1)
            {
                buf[5]++;
                //MPI::COMM_WORLD.Send(&n,1,mpiInt,top_rank,tag);
                r1[5] = mpiIRecv(&rec[1],1,mpiInt,top_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< top_rank  << endl;
            }
            if(bottom_rank!= 1)
            {
                buf[6]++;
                //MPI::COMM_WORLD.Send(&n,1,mpiInt,bottom_rank,tag);
                r1[6] = mpiIRecv(&rec[2],1,mpiInt,bottom_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< bottom_rank  << endl;
            }
            if(right_rank!= 1)
            {
                buf[7]++;
                //MPI::COMM_WORLD.Send(&n,1,mpiInt,right_rank,tag);
                r1[7] = mpiIRecv(&rec[3],1,mpiInt,right_rank,tag);
                //cout << "This is "<< current_rank << " received "<< rec << " from "<< right_rank <<endl;
            }
        }
        for(int i =0; i< 8; i++)
        {
            if(buf[i] !=0)
            {
                r1[i].Wait();
            }
        }
//cout << rec[0] <<" "<< rec[1]<<" " << rec[2]<< " " << rec[3] << endl;
        /*
        else{
        //cout << " current rank is"<< current_rank << endl;
        	if(left_rank != 1){
        	//MPI::COMM_WORLD.Send(&n,1,mpiInt,left_rank,tag);
        	r1 = mpiIRecv(&rec,1,mpiInt,left_rank,tag);
        	r1.Wait();
        	cout << "This is "<< current_rank << " received "<< rec << " from "<< left_rank <<endl ;
        	}
        	if(top_rank != 1){
        	//MPI::COMM_WORLD.Send(&n,1,mpiInt,top_rank,tag);
        	r1=mpiIRecv(&rec,1,mpiInt,top_rank,tag);
        	r1.Wait();
        	cout << "This is "<< current_rank << " received "<< rec << " from "<< top_rank  << endl;
        	}
        	if(bottom_rank != 1){
        	//MPI::COMM_WORLD.Send(&n,1,mpiInt,bottom_rank,tag);
        	r1=mpiIRecv(&rec,1,mpiInt,bottom_rank,tag);
        	r1.Wait();
        	cout << "This is "<< current_rank << " received "<< rec << " from "<< bottom_rank  << endl;
        	}
        	if(right_rank != 1){
        	//MPI::COMM_WORLD.Send(&n,1,mpiInt,right_rank,tag);
        	r1=mpiIRecv(&rec,1,mpiInt,right_rank,tag);
        	r1.Wait();
        	cout << "This is "<< current_rank << " received "<< rec << " from "<< right_rank <<endl;
        	}

        	if(left_rank!= 1){
        	r1=mpiISend(&n,1,mpiInt,left_rank,tag);
        	//MPI::COMM_WORLD.Recv(&rec,1,mpiInt,left_rank,tag);
        	//cout << "This is "<< current_rank << " received "<< rec << " from "<< left_rank <<endl ;
        	}
        	if(top_rank!= 1){
        	r1=mpiISend(&n,1,mpiInt,top_rank,tag);
        	//MPI::COMM_WORLD.Recv(&rec,1,mpiInt,top_rank,tag);
        	//cout << "This is "<< current_rank << " received "<< rec << " from "<< top_rank  << endl;
        	}
        	if(bottom_rank!= 1){
        	r1=mpiISend(&n,1,mpiInt,bottom_rank,tag);
        	//MPI::COMM_WORLD.Recv(&rec,1,mpiInt,bottom_rank,tag);
        	//cout << "This is "<< current_rank << " received "<< rec << " from "<< bottom_rank  << endl;
        	}
        	if(right_rank!= 1){
        	r1=mpiISend(&n,1,mpiInt,right_rank,tag);
        	//MPI::COMM_WORLD.Recv(&rec,1,mpiInt,right_rank,tag);
        	//cout << "This is "<< current_rank << " received "<< rec << " from "<< right_rank <<endl;
        	}

        }

        */

    }
//cout << "final "<< rec << endl;
}


void domain::Send_message(dataType *A,int size,int dest, int tag)
{
    MPI::COMM_WORLD.Send(A, size,mpiDType,dest,tag);
}


void domain::Recv_message(dataType *A,int max_size, int source, int tag)
{
    MPI::Status status;
    MPI::COMM_WORLD.Recv(A,max_size,mpiDType,source,tag,status);
}


void domain::Set_params(int a, int b, int c, int d)
{
// sets up initail parameters
    width = a;
    height = b;
    horz_proc = c;
    vert_proc = d;
    return;
}


void domain::Set_ranks(int rank,int size)
{
// set ranks of surroundings
    if(size != horz_proc*vert_proc+1)
    {
        //handling wrong number of process specification
        cout << "Wrong number of processes given, it should be "<< horz_proc*vert_proc+1 << endl;
        return;
    }
    current_rank = rank; //setting rank of the current domain

    if(current_rank == 0)
    {
        //master process
        left_rank =0;
        right_rank =0;
        top_rank =0;
        bottom_rank =0;
    }
    //else if(current_rank == 1){
    // boundary process
    // 	left_rank =1; right_rank =1; top_rank =1; bottom_rank =1;
    //}
    else
    {
        // remaining other processes
        int proc_row,proc_column;
        proc_row = (rank-1)%vert_proc; // process row
        proc_column = floor((rank-1)/vert_proc); // process column
        top_rank = current_rank - 1;
        right_rank = current_rank +vert_proc;
        bottom_rank = current_rank + 1;
        left_rank = current_rank - vert_proc;
        if(proc_row == 0)
        {
            top_rank =-1; // wall boundary
        }
        if(proc_row == vert_proc-1)
        {
            bottom_rank = -1; // wall boundary
        }
        if(proc_column == 0)
        {
            left_rank = -2;	// inflow
        }
        if(proc_column == horz_proc -1)
        {
            right_rank = -3; // outflow
        }

    }
    swap(&top_rank,&bottom_rank);
    //cout << "currnt process is " << current_rank << " with "<< left_rank << " " << top_rank << " " << right_rank << " " << bottom_rank << " " <<endl;
    return;
}



void domain::Set_indexes()
{
// sets the domain grid indexes
    if(current_rank == 0)
    {
        //master process
        grid_bottom.x =0;
        grid_bottom.y =0;
        grid_top.x =0;
        grid_top.y=0;
    }
    //else if(current_rank == 1){
    // boundary process
    // 	grid_bottom.x =1; grid_bottom.y =1; grid_top.x =1; grid_top.y=1;
    //}
    else
    {
        // remaining other processes
        int proc_row,proc_column,row_points,column_points;
        proc_row = (current_rank-1)%vert_proc; // process row
        proc_column = floor((current_rank-1)/vert_proc); // process column
        row_points = floor(height/vert_proc);
        column_points = floor(width/horz_proc);
        grid_bottom.y = proc_row*row_points;
        grid_top.y = grid_bottom.y + row_points - 1;
        grid_bottom.x = proc_column *column_points;
        grid_top.x = grid_bottom.x + column_points -1;
        if(proc_row == vert_proc-1)
        {
            grid_top.y = height - 1;
        }
        if(proc_column == horz_proc -1)
        {
            grid_top.x = width -1;
        }
//cout << "currnt process is " << current_rank << " with "<< grid_bottom.x << " " << grid_bottom.y << " "<< grid_top.x << " "<< grid_top.y << " " << endl;
    }
    return;
}


class grid
{
// load grid nodes
};



int main(int argc, char *argv[])
{
    domain params; // domain or process params
    params.Set_params(97,49,1,1);
    int total_proc = params.horz_proc * params.vert_proc + 1; // one to handle boundary and other is parent process
    int rank,size;
    MPI::Init(); // initializing MPI
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    params.Load_grid();
    params.Set_ranks(rank,size);
    params.Set_indexes();
    for(int i =2; i < 3; i++)
    {
//params.Send_to_surr(pow(rank,i));
    }
//cout << MPI::TAG_UB <<endl;
    MPI::Finalize(); // end MPI
    return 0;
}
