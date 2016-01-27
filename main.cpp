/**
  *
  * parallel flow solver using MPI
  * @authr R D Teja - rdarmateja@gmail.com
  */
#include <mpi.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <sstream>

/* patch to make to_string work */
namespace patch
{
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}

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
    dataType R;//gas constant
    dataType r;//gaama
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
    ///Main class for the process, it sets all the required parameters to work with
    /// domain parameters, grid nodes, num of processes in horiz and vertical
public:
    int width; // horizontal grid points
    int height; // vertical grid points
    int horz_proc;// horizontal processes (cols)
    int vert_proc;// vertical processes (rows)
    ///so total processes will be (h*v + 2)
    int left_rank,top_rank,right_rank,bottom_rank,current_rank;//ranks of surrounding processes computed by Set_ranks
    grid_index grid_bottom,grid_top; //stores the index of the grid points that the process is operating on, computed by set_indexes
    /// grid data X
    std::vector <dataType> X;
    std::vector <dataType> Y;
    /// flow prop
    flow_variables in_flow;
    /// prototype declarations
    void Set_params(int a, int b, int c, int d);// set the parameters
    void Set_ranks(int rank,int size); // sets surrounding processes
    void Set_indexes(); // sets grid indexres
    void Send_message(dataType *A,int size,int dest, int tag);
    void Recv_message(dataType *A,int max_size, int source, int tag);
    void Send_to_surr(int n);
    void Send_to_surr2();
    void Load_grid(std::string str);
    void Run_vanleer();//van-leer scheme, 2d
    void Write_coord(); // write coordinates to text file
    /*
    order of running
    1. Set_params()
    2. Load_grid()
    3. Set_ranks()
    4. Set_indexes()
    5. Run_vanleer()
    */
};

void domain::Run_vanleer()
{
    /// vanleer scheme
    int imx = width, jmx = height;

    /******************storing coordinates and computing areas, normals, volume start*********************/
    dataType x_cord[imx][jmx],y_cord[imx][jmx];//x and y coordinates of grid
    dataType Iar[imx][jmx],Jar[imx][jmx];//i-face and j-face areas
    dataType In[imx][jmx][2],Jn[imx][jmx][2],l;// i,j face normals
    dataType Vol[imx][jmx],lo,up;// volume of cells
    /// Storing the coordinates
    int buf =0;
    for(int j=0 ; j<jmx ; j++)
    {
        for(int i=0 ; i<imx ; i++)
        {
            x_cord[i][j] = X[buf];
            y_cord[i][j] = Y[buf];
            buf++;
        }
    }
    /// Calculating i face areas
    for(int j=0; j<jmx-1; j++)
    {
        for(int i=0; i<imx; i++)
        {
            Iar[i][j+1] = sqrt(pow((x_cord[i][j+1] - x_cord[i][j]),2)+ pow((y_cord[i][j+1] - y_cord[i][j]),2));
        }
    }

    /// Calculating j face areas
    for(int j=0; j<jmx; j++)
    {
        for(int i=0; i<imx-1; i++)
        {
            Jar[i+1][j] = sqrt(pow((x_cord[i+1][j] - x_cord[i][j]),2)+ pow((y_cord[i+1][j] - y_cord[i][j]),2));
        }
    }
    /// Calculating i face normals
    for(int j=0; j<jmx-1; j++)
    {
        for(int i=0; i<imx; i++)
        {
            l = sqrt(pow((x_cord[i][j+1] - x_cord[i][j]),2)+ pow((y_cord[i][j+1] - y_cord[i][j]),2));
            In[i][j+1][0] = (y_cord[i][j+1] - y_cord[i][j])/l; // x-component of i-face normal
            In[i][j+1][1] = (-x_cord[i][j+1] + x_cord[i][j])/l; // y-component of i-face normal
        }
    }

    /// Calculating j face normals
    for(int j=0; j<jmx; j++)
    {
        for(int i=0; i<imx-1; i++)
        {
            l = sqrt(pow((x_cord[i+1][j] - x_cord[i][j]),2) + pow((y_cord[i+1][j] - y_cord[i][j]),2))	;
            Jn[i+1][j][0] = (-y_cord[i+1][j] + y_cord[i][j])/l; // x-component of j-face normal
            Jn[i+1][j][1] = (x_cord[i+1][j] - x_cord[i][j])/l; // y-component of j-face normal
        }
    }
    /// Calculating cell volumes
    for(int j=0; j<jmx-1; j++)
    {
        for(int i=0; i<imx-1; i++)
        {
            lo = 0.5*fabs((x_cord[i+1][j] - x_cord[i][j])*(y_cord[i+1][j+1] - y_cord[i+1][j]) - (x_cord[i+1][j+1] - x_cord[i+1][j])*(y_cord[i+1][j] - y_cord[i][j]));
            up = 0.5*fabs((x_cord[i][j+1] - x_cord[i][j])*(y_cord[i+1][j+1] - y_cord[i][j+1]) - (x_cord[i+1][j+1] - x_cord[i][j+1])*(y_cord[i][j+1] - y_cord[i][j]));
            Vol[i+1][j+1] = lo + up; // Volume of each cell
        }
    }

    /******************storing coordinates and computing areas,normal,volumes end*********************/

    /// inflow conditions
    dataType pInf = in_flow.press,dInf = in_flow.dens, tInf = in_flow.temp;
    dataType r = in_flow.r, R = in_flow.R;
    dataType mInf = in_flow.mach,uInf = in_flow.u,vInf = in_flow.v,wInf = in_flow.w;

}


void domain::Write_coord()
{
    /// getting coordinates
    int imx = width, jmx = height;
    dataType x_cord[imx][jmx],y_cord[imx][jmx];//x and y coordinates of grid
    int buf =0;
    for(int j=0 ; j<jmx ; j++)
    {
        for(int i=0 ; i<imx ; i++)
        {
            x_cord[i][j] = X[buf];
            y_cord[i][j] = Y[buf];
            buf++;
        }
    }
    /// writes process coordinates to a text file
    string str = "txt/process_coord_" + patch::to_string(current_rank) + ".txt";
    char *c = &str[0u];
    std::fstream coord (c,fstream::out);
    if(coord.is_open())
    {
        for(int j = grid_bottom.y; j <= grid_top.y; j++)
        {
            for(int i = grid_bottom.x; i <= grid_top.x; i++)
            {
                coord << x_cord[i][j] << " " << y_cord[i][j] << "\n";
            }
        }
        coord.close();
    }

}


void domain::Load_grid(std::string str)
{
// to load grid into an array
    char *c = &str[0u];
    std::ifstream grid;
    grid.open(c);
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
    grid.close();
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
    /*
    else if(current_rank == 1){
     boundary process
     	left_rank =1; right_rank =1; top_rank =1; bottom_rank =1;
    }
    */
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
    swap(&top_rank,&bottom_rank); // just a small manipulation
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
    /*
    else if(current_rank == 1){
     boundary process
     	grid_bottom.x =1; grid_bottom.y =1; grid_top.x =1; grid_top.y=1;
    }
    */
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

    int rank,size;
    MPI::Init(); // initializing MPI
    size = MPI::COMM_WORLD.Get_size();
    rank = MPI::COMM_WORLD.Get_rank();
    if(rank != 0)
    {
        domain params; // domain or process params
        params.Set_params(0,0,1,1);
        int total_proc = params.horz_proc * params.vert_proc + 1; // one to handle boundary and other is parent process
        params.Load_grid("bumpgrid.txt");// laod grid
        params.Set_ranks(rank,size); // set ranks of surrounding processes
        params.Set_indexes(); // compute indexes of the process
        /*set inflow conditions*/
        params.in_flow.press = 1013.26;
        params.in_flow.temp = 300;
        params.in_flow.R = 287;
        params.in_flow.r = 1.4;
        params.in_flow.mach = 0.4;
        params.in_flow.dens = params.in_flow.press/(params.in_flow.R*params.in_flow.temp);
        params.in_flow.u = params.in_flow.mach*sqrt(params.in_flow.r*params.in_flow.press/params.in_flow.dens);
        params.in_flow.v =0;
        params.in_flow.w =0;
        //cout << params.in_flow.dens << endl;
        params.Run_vanleer();
        //params.Write_coord();
        for(int i =2; i < 10; i++)
        {
          //  params.Send_to_surr(pow(rank,i));
        }
        //cout << MPI::TAG_UB <<endl;
    }
    MPI::Finalize(); // end MPI
    return 0;
}
