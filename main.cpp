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

#define dataType long double
#define mpiDType MPI::LONG_DOUBLE
#define mpiInt MPI::INT
#define mpiSend MPI::COMM_WORLD.Send
#define mpiRecv MPI::COMM_WORLD.Recv
#define mpiISend MPI::COMM_WORLD.Irecv
#define mpiIRecv MPI::COMM_WORLD.Isend


namespace patch
{
/* patch to make to_string work */
template < typename T > std::string to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}
}

/*supports for vanleer scheme*/
namespace vanleer
{
/*supports for vanleer scheme*/

/// Defining alpha +
double alp(double M)
{
    if (M > 0)
        return(1);
    return(0);
}

int sign(double x)
{
   return (x>0)- (x<0);
}

}


using namespace std;
/******************** supports start************************/
struct  flow_variables
{
    /* data */
    int iterations;
    dataType cfl;
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
    6. Write_coord()
    */
};

void domain::Run_vanleer()
{
    /// vanleer scheme
    int imx = width, jmx = height; // imax and jmax of entire flow domain

    /******************storing coordinates and computing areas, normals, volume start*********************/
    dataType x_cord[imx][jmx],y_cord[imx][jmx];//x and y coordinates of grid
    dataType Iar[imx][jmx-1],Jar[imx-1][jmx];//i-face and j-face areas
    dataType In[imx][jmx-1][2],Jn[imx-1][jmx][2],l;// i,j face normals
    dataType Vol[imx][jmx],lo,up;// volume of cells
    /// Storing the coordinates of entire domain
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
    /// getting coordinates of process domain
    /*
    imx = grid_top.x-grid_bottom.x+1;
    jmx = grid_top.y-grid_bottom.y+1;
    dataType x_cord2[imx][jmx],y_cord2[imx][jmx];
    dataType Iar[imx][jmx],Jar[imx][jmx];//i-face and j-face areas
    dataType In[imx][jmx][2],Jn[imx][jmx][2],l;// i,j face normals
    dataType Vol[imx][jmx],lo,up;// volume of cells

    for(int j = grid_bottom.y; j <= grid_top.y; j++)
    {
        for(int i = grid_bottom.x; i <= grid_top.x; i++)
        {
            //coord << x_cord[i][j] << " " << y_cord[i][j] << "\n";
            x_cord2[i-grid_bottom.x][j-grid_bottom.y] = x_cord[i][j];
            y_cord2[i-grid_bottom.x][j-grid_bottom.y] = y_cord[i][j];
        }
    }
    */
    /// Calculating i face areas
    for(int j=0; j<jmx-1; j++)
    {
        for(int i=0; i<imx; i++)
        {
            Iar[i][j] = sqrt(pow((x_cord[i][j+1] - x_cord[i][j]),2)+ pow((y_cord[i][j+1] - y_cord[i][j]),2));
        }
    }

    /// Calculating j face areas
    for(int j=0; j<jmx; j++)
    {
        for(int i=0; i<imx-1; i++)
        {
            Jar[i][j] = sqrt(pow((x_cord[i+1][j] - x_cord[i][j]),2)+ pow((y_cord[i+1][j] - y_cord[i][j]),2));
        }
    }
    /// Calculating i face normals
    for(int j=0; j<jmx-1; j++)
    {
        for(int i=0; i<imx; i++)
        {
            l = sqrt(pow((x_cord[i][j+1] - x_cord[i][j]),2)+ pow((y_cord[i][j+1] - y_cord[i][j]),2));
            In[i][j][0] = (y_cord[i][j+1] - y_cord[i][j])/l; // x-component of i-face normal
            In[i][j][1] = (-x_cord[i][j+1] + x_cord[i][j])/l; // y-component of i-face normal
        }
    }

    /// Calculating j face normals
    for(int j=0; j<jmx; j++)
    {
        for(int i=0; i<imx-1; i++)
        {
            l = sqrt(pow((x_cord[i+1][j] - x_cord[i][j]),2) + pow((y_cord[i+1][j] - y_cord[i][j]),2))	;
            Jn[i][j][0] = (-y_cord[i+1][j] + y_cord[i][j])/l; // x-component of j-face normal
            Jn[i][j][1] = (x_cord[i+1][j] - x_cord[i][j])/l; // y-component of j-face normal
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
    /// primitive and Characteristic variables decleration
    int imax= grid_top.x-grid_bottom.x+1, jmax = grid_top.y-grid_bottom.y+1; // imax and jmax of current domain
    //cout << imax << " "<<jmax << endl;
    dataType q[imax+1][jmax+1][4]; // primitive, size includes ghost cells for boundary,in-outflow, interface
    dataType U[imax+1][jmax+1][4]; // characteristic
    /// initializing to inflow conditions
    for(int i =0; i < imax+1; i++)
    {
        for(int j =0; j< jmax+1; j++)
        {
            q[i][j][0] = pInf;
            q[i][j][1] = uInf;
            q[i][j][2] = vInf;
            q[i][j][3] = dInf;
            U[i][j][0] = dInf;
            U[i][j][1] = dInf*uInf;
            U[i][j][2] = dInf*vInf;
            U[i][j][3] = pInf/(r-1) + 0.5*dInf*(uInf*uInf + vInf*vInf);
        }
    }
    /****************************apply b.c. outsite main loop - starts***********************************/
    /// left face - can be inflow or interface(no need to impose any conditions)
    if(left_rank == -2)
    {
        /// subsonic
        if(mInf < 1)
        {
            for(int j =0; j< jmax+1; j++)
            {
                q[0][j][0] = q[1][j][0];// pressure is imposed from right
            }
        }
        /// supersonic - implicitly satisfied as everything is initialized to infi values
    }
    /// top face - can be boundary or interface
    if(top_rank == -1)
    {
        for(int i =0; i < imax+1; i++)
        {
            q[i][jmax][0] =  q[i][jmax-1][0];// pressure
            q[i][jmax][3] =  q[i][jmax-1][3];// density
        }

    }

    /// right face - can be outflow or interface
    if(right_rank == -3)
    {
        /*
        if(mInf < 1)
        {
            for(int j=0; j<jmax+1; j++)
            {
            }

        }
        */
    }
    /// bottom face - can be boundary or interface
    if(bottom_rank == -1)
    {
        for(int i =0; i < imax+1; i++)
        {
            q[i][0][0] =  q[i][1][0];// pressure
            q[i][0][3] =  q[i][1][3];// density
        }
    }
    /****************************apply b.c. outsite main loop - ends***********************************/

    /**************************************main iteration loop - starts ***************************/
    /// switch variable
    dataType mp1[imax][jmax-1],mm1[imax][jmax-1],mp2[imax-1][jmax],mm2[imax-1][jmax];//mach + ,-
    dataType ap1[imax][jmax-1],am1[imax][jmax-1],ap2[imax-1][jmax],am2[imax-1][jmax];//alpha + ,-
    dataType bp1[imax][jmax-1],bm1[imax][jmax-1],bp2[imax-1][jmax],bm2[imax-1][jmax];//beta + ,-
    dataType cp1[imax][jmax-1],cm1[imax][jmax-1],cp2[imax-1][jmax],cm2[imax-1][jmax];//cvl + ,-
    dataType dp1[imax][jmax-1],dm1[imax][jmax-1],dp2[imax-1][jmax],dm2[imax-1][jmax];//d + ,-
    dataType a[imax+1][jmax+1];// sound speed
    /// fluxes
    dataType F[imax][jmax-1][4],G[imax-1][jmax][4];
    /// residuals
    dataType Res[imax+1][jmax+1][4];
    /// time step
    dataType cfl = in_flow.cfl;
    dataType vnp1[imax][jmax-1],vnm1[imax][jmax-1],vnp2[imax-1][jmax],vnm2[imax-1][jmax];
    dataType rp1[imax][jmax-1],rm1[imax][jmax-1],rp2[imax-1][jmax],rm2[imax-1][jmax];
    dataType tvl[imax+1][jmax+1];
    /// message passing variables - sends and receives  concatinated primitive variables
    int left_len = (jmax+1)*4, top_len = (imax+1)*4;
    dataType Left_send[left_len],Right_send[left_len],Left_recv[left_len],Right_recv[left_len];
    dataType Top_send[top_len],Bottom_send[top_len],Top_recv[top_len],Bottom_recv[top_len];
    int counter;
    int iters = 0;
    /// resnorm
    dataType resNorm[in_flow.iterations];
    while(++iters <= in_flow.iterations)
    {
        resNorm[iters-1] =0;
        ///cout << "teja" << endl;
        /// sound speed
        for(int i =0; i< imax+1; i++)
        {
            for(int j =0; j< jmax+1; j++)
            {
                ///a(1:imax+1,1:jmax+1) = (sqrt(rair*q(1:imax+1,1:jmax+1,1)./q(1:imax+1,1:jmax+1,4))); % sound speeds
                a[i][j] = sqrt(r*q[i][j][0]/q[i][j][3]);
        }
        }
        //cout << dInf << endl;
        /// mach, alpha , beta, c, d for i faces
        int offset_i, offset_j;
        for(int j =0; j< jmax-1; j++)
        {
            for(int i =0; i< imax; i++ )
            {
                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                mp1[i][j] = (q[i][j+1][1]*In[offset_i][offset_j][0] + q[i][j+1][2]*In[offset_i][offset_j][1])/(0.5*(a[i][j+1]+a[i+1][j+1]));
                mm1[i][j] = (q[i+1][j+1][1]*In[offset_i][offset_j][0] + q[i+1][j+1][2]*In[offset_i][offset_j][1])/(0.5*(a[i][j+1]+a[i+1][j+1]));
                bp1[i][j] = std::min(0,vanleer::sign(std::abs(mp1[i][j])-1));
                bm1[i][j] = std::min(0,vanleer::sign(std::abs(mm1[i][j])-1));
                ap1[i][j] = 0.5*(1+vanleer::sign(mp1[i][j]));
                am1[i][j] = 0.5*(1-vanleer::sign(mm1[i][j]));
                ///cp1 =ap1.*(1+bp1).*mp1 - 0.25*bp1.*((1+mp1).^2);
                ///cm1 =am1.*(1+bm1).*mm1 + 0.25*bm1.*((1-mm1).^2);
                cp1[i][j] = ap1[i][j]*(1+bp1[i][j])*mp1[i][j] - 0.25*bp1[i][j]*(1+mp1[i][j])*(1+mp1[i][j]);
                cm1[i][j] = am1[i][j]*(1+bm1[i][j])*mm1[i][j] + 0.25*bm1[i][j]*(1-mm1[i][j])*(1-mm1[i][j]);
                /// dp1 =ap1.*(1+bp1)-0.25*bp1.*((1+mp1).^2).*(2-mp1);
                ///dm1 =am1.*(1+bm1)-0.25*bm1.*((1-mm1).^2).*(2+mm1);
                dp1[i][j] = ap1[i][j]*(1+bp1[i][j]) - 0.25*bp1[i][j]*(1+mp1[i][j])*(1+mp1[i][j])*(2-mp1[i][j]);
                dm1[i][j] = am1[i][j]*(1+bm1[i][j]) - 0.25*bm1[i][j]*(1-mm1[i][j])*(1-mm1[i][j])*(2+mm1[i][j]);
            }
        }

        /// mach, alpha , beta, c, d for j faces
        for(int j =0; j< jmax; j++)
        {
            for(int i =0; i< imax-1; i++ )
            {
                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                mp2[i][j] = (q[i+1][j][1]*Jn[offset_i][offset_j][0] + q[i+1][j][2]*Jn[offset_i][offset_j][1])/(0.5*(a[i+1][j]+a[i+1][j+1]));
                mm2[i][j] = (q[i+1][j+1][1]*Jn[offset_i][offset_j][0] + q[i+1][j+1][2]*Jn[offset_i][offset_j][1])/(0.5*(a[i+1][j]+a[i+1][j+1]));
                bp2[i][j] = std::min(0,vanleer::sign(std::abs(mp2[i][j])-1));
                bm2[i][j] = std::min(0,vanleer::sign(std::abs(mm2[i][j])-1));
                ap2[i][j] = 0.5*(1+vanleer::sign(mp2[i][j]));
                am2[i][j] = 0.5*(1-vanleer::sign(mm2[i][j]));
                ///cp1 =ap1.*(1+bp1).*mp1 - 0.25*bp1.*((1+mp1).^2);
                ///cm1 =am1.*(1+bm1).*mm1 + 0.25*bm1.*((1-mm1).^2);
                cp2[i][j] = ap2[i][j]*(1+bp2[i][j])*mp2[i][j] - 0.25*bp2[i][j]*(1+mp2[i][j])*(1+mp2[i][j]);
                cm2[i][j] = am2[i][j]*(1+bm2[i][j])*mm2[i][j] + 0.25*bm2[i][j]*(1-mm2[i][j])*(1-mm2[i][j]);
                /// dp1 =ap1.*(1+bp1)-0.25*bp1.*((1+mp1).^2).*(2-mp1);
                ///dm1 =am1.*(1+bm1)-0.25*bm1.*((1-mm1).^2).*(2+mm1);
                dp2[i][j] = ap2[i][j]*(1+bp2[i][j]) - 0.25*bp2[i][j]*(1+mp2[i][j])*(1+mp2[i][j])*(2-mp2[i][j]);
                dm2[i][j] = am2[i][j]*(1+bm2[i][j]) - 0.25*bm2[i][j]*(1-mm2[i][j])*(1-mm2[i][j])*(2+mm2[i][j]);

            }
        }

        /// F fluxes
        for(int j=0; j<jmax-1; j++)
        {
            for(int i =0; i< imax; i++)
            {

                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                /*
                    fFlux1(1:imax,1:jmax-1,1) = (q(1:imax,2:jmax,4).*a(1:imax,2:jmax).*cp1 + q(2:imax+1,2:jmax,4).*a(2:imax+1,2:jmax).*cm1).*xareamag;
                fFlux1(1:imax,1:jmax-1,2) = (U(1:imax,2:jmax,2).*a(1:imax,2:jmax).*cp1 + U(2:imax+1,2:jmax,2).*a(2:imax+1,2:jmax).*cm1 + (dp1.*q(1:imax,2:jmax,1)+dm1.*q(2:imax+1,2:jmax,1)).*xarean(:,:,1)).*xareamag;
                fFlux1(1:imax,1:jmax-1,3) = (U(1:imax,2:jmax,3).*a(1:imax,2:jmax).*cp1 + U(2:imax+1,2:jmax,3).*a(2:imax+1,2:jmax).*cm1 + (dp1.*q(1:imax,2:jmax,1)+dm1.*q(2:imax+1,2:jmax,1)).*xarean(:,:,2)).*xareamag;
                fFlux1(1:imax,1:jmax-1,4) = (q(1:imax,2:jmax,4).*a(1:imax,2:jmax).*cp1.*(3.5*q(1:imax,2:jmax,1)./q(1:imax,2:jmax,4)+0.5*(q(1:imax,2:jmax,2).^2+q(1:imax,2:jmax,3).^2))...
                 + q(2:imax+1,2:jmax,4).*a(2:imax+1,2:jmax).*cm1.*(3.5*q(2:imax+1,2:jmax,1)./q(2:imax+1,2:jmax,4)+0.5*(q(2:imax+1,2:jmax,2).^2+q(2:imax+1,2:jmax,3).^2))).*xareamag;
                */
                F[i][j][0] = (q[i][j+1][3]*a[i][j+1]*cp1[i][j] + q[i+1][j+1][3]*a[i+1][j+1]*cm1[i][j])*Iar[offset_i][offset_j];
                F[i][j][1] = (q[i][j+1][3]*q[i][j+1][1]*a[i][j+1]*cp1[i][j] + q[i+1][j+1][3]*q[i+1][j+1][1]*a[i+1][j+1]*cm1[i][j] +
                (dp1[i][j]*q[i][j+1][0] + dm1[i][j]*q[i+1][j+1][0])*In[offset_i][offset_j][0])*Iar[offset_i][offset_j];
                F[i][j][2] = (q[i][j+1][3]*q[i][j+1][2]*a[i][j+1]*cp1[i][j] + q[i+1][j+1][3]*q[i+1][j+1][2]*a[i+1][j+1]*cm1[i][j] +
                (dp1[i][j]*q[i][j+1][0] + dm1[i][j]*q[i+1][j+1][0])*In[offset_i][offset_j][1])*Iar[offset_i][offset_j];
                F[i][j][3] = (q[i][j+1][3]*a[i][j+1]*cp1[i][j]*(3.5*q[i][j+1][0]/q[i][j+1][3] + 0.5*(q[i][j+1][1]*q[i][j+1][1] + q[i][j+1][2]*q[i][j+1][2])) +
                q[i+1][j+1][3]*a[i+1][j+1]*cm1[i][j]*(3.5*q[i+1][j+1][0]/q[i+1][j+1][3] + 0.5*(q[i+1][j+1][1]*q[i+1][j+1][1] + q[i+1][j+1][2]*q[i+1][j+1][2])))*Iar[offset_i][offset_j];
            }

        }
        /// G fluxes
        for(int j=0; j<jmax; j++)
        {
            for(int i =0; i< imax-1; i++)
            {

                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                /*
                    gFlux1(1:imax-1,1:jmax,1) = (q(2:imax,1:jmax,4).*a(2:imax,1:jmax).*cp2 + q(2:imax,2:jmax+1,4).*a(2:imax,2:jmax+1).*cm2).*yareamag;
                gFlux1(1:imax-1,1:jmax,2) = (U(2:imax,1:jmax,2).*a(2:imax,1:jmax).*cp2 + U(2:imax,2:jmax+1,2).*a(2:imax,2:jmax+1).*cm2 + (dp2.*q(2:imax,1:jmax,1)+dm2.*q(2:imax,2:jmax+1,1)).*yarean(:,:,1)).*yareamag;
                gFlux1(1:imax-1,1:jmax,3) = (U(2:imax,1:jmax,3).*a(2:imax,1:jmax).*cp2 + U(2:imax,2:jmax+1,3).*a(2:imax,2:jmax+1).*cm2 + (dp2.*q(2:imax,1:jmax,1)+dm2.*q(2:imax,2:jmax+1,1)).*yarean(:,:,2)).*yareamag;
                gFlux1(1:imax-1,1:jmax,4) = (q(2:imax,1:jmax,4).*a(2:imax,1:jmax).*cp2.*(3.5*q(2:imax,1:jmax,1)./q(2:imax,1:jmax,4)+0.5*(q(2:imax,1:jmax,2).^2+q(2:imax,1:jmax,3).^2))...
                + q(2:imax,2:jmax+1,4).*a(2:imax,2:jmax+1).*cm2.*(3.5*q(2:imax,2:jmax+1,1)./q(2:imax,2:jmax+1,4)+0.5*(q(2:imax,2:jmax+1,2).^2+q(2:imax,2:jmax+1,3).^2))).*yareamag;

                */
                G[i][j][0] = (q[i+1][j][3]*a[i+1][j]*cp2[i][j] + q[i+1][j+1][3]*a[i+1][j+1]*cm2[i][j])*Jar[offset_i][offset_j];
                G[i][j][1] = (q[i+1][j][3]*q[i+1][j][1]*a[i+1][j]*cp2[i][j] + q[i+1][j+1][3]*q[i+1][j+1][1]*a[i+1][j+1]*cm2[i][j] +
                (dp2[i][j]*q[i+1][j][0] + dm2[i][j]*q[i+1][j+1][0])*Jn[offset_i][offset_j][0])*Jar[offset_i][offset_j];
                G[i][j][2] = (q[i+1][j][3]*q[i+1][j][2]*a[i+1][j]*cp2[i][j] + q[i+1][j+1][3]*q[i+1][j+1][2]*a[i+1][j+1]*cm2[i][j] +
                (dp2[i][j]*q[i+1][j][0] + dm2[i][j]*q[i+1][j+1][0])*Jn[offset_i][offset_j][1])*Jar[offset_i][offset_j];
                G[i][j][3] = (q[i+1][j][3]*a[i+1][j]*cp2[i][j]*(3.5*q[i+1][j][0]/q[i+1][j][3] + 0.5*(q[i+1][j][1]*q[i+1][j][1] + q[i+1][j][2]*q[i+1][j][2])) +
                q[i+1][j+1][3]*a[i+1][j+1]*cm2[i][j]*(3.5*q[i+1][j+1][0]/q[i+1][j+1][3] + 0.5*(q[i+1][j+1][1]*q[i+1][j+1][1] + q[i+1][j+1][2]*q[i+1][j+1][2])))*Jar[offset_i][offset_j];
            }

        }
        /********************** wall no penetration boundary conditions on G flux - starts********************/
        /// top wall
        if(top_rank == -1)
        {
            for(int i = 0; i <  imax-1; i++)
            {
                offset_i = i+ grid_bottom.x;
                offset_j = grid_bottom.y;
                G[i][jmax-1][0] = 0;
                G[i][jmax-1][3] = 0;
                /*
                    gFlux1(1:imax-1,jmax,2)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,1);
                    gFlux1(1:imax-1,jmax,3)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,2);
                */
                G[i][jmax-1][1] = (dp2[i][jmax-1]*q[i+1][jmax][0] + dm2[i][jmax-1]*q[i+1][jmax][0])*Jn[offset_i][offset_j+jmax-1][0]*Jar[offset_i][offset_j+jmax-1];
                G[i][jmax-1][2] = (dp2[i][jmax-1]*q[i+1][jmax][0] + dm2[i][jmax-1]*q[i+1][jmax][0])*Jn[offset_i][offset_j+jmax-1][1]*Jar[offset_i][offset_j+jmax-1];
            }
        }
        if(bottom_rank == -1)
        {
            for(int i = 0; i <  imax-1; i++)
            {
                offset_i = i+ grid_bottom.x;
                offset_j = grid_bottom.y;
                G[i][0][0] = 0;
                G[i][0][3] = 0;
                /*
                    gFlux1(1:imax-1,jmax,2)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,1);
                    gFlux1(1:imax-1,jmax,3)= (dp2(1:imax-1,jmax).*q(2:imax,jmax+1,1)+dm2(1:imax-1,jmax).*q(2:imax,jmax+1,1)).*yarea(1:imax-1,jmax,2);
                */
                G[i][0][1] = (dp2[i][0]*q[i+1][1][0] + dm2[i][0]*q[i+1][1][0])*Jn[offset_i][offset_j][0]*Jar[offset_i][offset_j];
                G[i][0][2] = (dp2[i][0]*q[i+1][1][0] + dm2[i][0]*q[i+1][1][0])*Jn[offset_i][offset_j][1]*Jar[offset_i][offset_j];
            }

        }
        /********************** wall boundary conditions on G flux -  ends********************/

        /***********************computing residuals -start ********************/
        for(int j = 1; j < jmax; j++)
        {
            for(int i =1; i< imax ; i++)
            {
                Res[i][j][0] = F[i][j-1][0] - F[i-1][j-1][0] + G[i-1][j][0] - G[i-1][j-1][0];
                Res[i][j][1] = F[i][j-1][1] - F[i-1][j-1][1] + G[i-1][j][1] - G[i-1][j-1][1];
                Res[i][j][2] = F[i][j-1][2] - F[i-1][j-1][2] + G[i-1][j][2] - G[i-1][j-1][2];
                Res[i][j][3] = F[i][j-1][3] - F[i-1][j-1][3] + G[i-1][j][3] - G[i-1][j-1][3];
                resNorm[iters-1] += Res[i][j][0]*Res[i][j][0] + Res[i][j][1]*Res[i][j][1] +Res[i][j][2]*Res[i][j][2]+Res[i][j][3]*Res[i][j][3];
            }
        }
        resNorm[iters-1] = sqrt(resNorm[iters-1]);
        /***********************computing residuals -stop ********************/

        /************************Time step calculation - start****************/
        /*
        vnp1(1:imax,1:jmax-1) = q(1:imax,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+ q(1:imax,2:jmax,3).*xarean(1:imax,1:jmax-1,2);
        vnm1(1:imax,1:jmax-1) = q(2:imax+1,2:jmax,2).*xarean(1:imax,1:jmax-1,1)+ q(2:imax+1,2:jmax,3).*xarean(1:imax,1:jmax-1,2);
        rp1 = abs(vnp1)+sqrt(rair*q(1:imax,2:jmax,1)./q(1:imax,2:jmax,4));
        rm1 = abs(vnm1)+sqrt(rair*q(2:imax+1,2:jmax,1)./q(2:imax+1,2:jmax,4));
        vnp2(1:imax-1,1:jmax) = q(2:imax,1:jmax,2).*yarean(1:imax-1,1:jmax,1)+ q(2:imax,1:jmax,3).*yarean(1:imax-1,1:jmax,2);
        rp2= abs(vnp2)+sqrt(rair*q(2:imax,1:jmax,1)./q(2:imax,1:jmax,4));
        vnm2(1:imax-1,1:jmax) = q(2:imax,2:jmax+1,2).*yarean(1:imax-1,1:jmax,1)+ q(2:imax,2:jmax+1,3).*yarean(1:imax-1,1:jmax,2);
        rm2 = abs(vnm2)+sqrt(rair*q(2:imax,2:jmax+1,1)./q(2:imax,2:jmax+1,4));
        */
        for(int j =0; j< jmax-1; j++)
        {
            for(int i =0; i< imax; i++)
            {
                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                vnp1[i][j] = q[i][j+1][1]*In[offset_i][offset_j][0] + q[i][j+1][2]*In[offset_i][offset_j][1];
                vnm1[i][j] = q[i+1][j+1][1]*In[offset_i][offset_j][0]+q[i+1][j+1][2]*In[offset_i][offset_j][1];
                rp1[i][j] = abs(vnp1[i][j])+sqrt(r*q[i][j+1][0]/q[i][j+1][3]);
                rm1[i][j] = abs(vnm1[i][j])+sqrt(r*q[i+1][j+1][0]/q[i+1][j+1][3]);
            }
        }

        for(int j =0; j< jmax; j++)
        {
            for(int i =0; i< imax-1; i++)
            {
                offset_i = i+grid_bottom.x;
                offset_j = j+grid_bottom.y;
                vnp2[i][j] = q[i+1][j][1]*Jn[offset_i][offset_j][0] + q[i+1][j][2]*Jn[offset_i][offset_j][1];
                vnm2[i][j] = q[i+1][j+1][1]*Jn[offset_i][offset_j][0]+q[i+1][j+1][2]*Jn[offset_i][offset_j][1];
                rp2[i][j] = abs(vnp2[i][j])+sqrt(r*q[i+1][j][0]/q[i+1][j][3]);
                rm2[i][j] = abs(vnm2[i][j])+sqrt(r*q[i+1][j+1][0]/q[i+1][j+1][3]);
            }
        }
         /// time step calculation

        /*


    tlv(2:imax,2:jmax) = (2*cfl)./(xareamag(1:imax-1,1:jmax-1).*rm1(1:imax-1,1:jmax-1) + xareamag(2:imax,1:jmax-1).*rp1(2:imax,1:jmax-1)
    +yareamag(1:imax-1,1:jmax-1).*rm2(1:imax-1,1:jmax-1)+yareamag(1:imax-1,2:jmax).*rp2(1:imax-1,2:jmax));

        */
        for(int j =1; j< jmax; j++ )
        {
            for(int i =1; i< imax; i++)
            {
                offset_i = i+grid_bottom.x-1;
                offset_j = j+grid_bottom.y-1;
                tvl[i][j] = (2*cfl)/(Iar[offset_i][offset_j]*rm1[i-1][j-1] + Iar[offset_i+1][offset_j]*rp1[i][j-1] +
                Jar[offset_i][offset_j]*rm2[i-1][j-1] + Jar[offset_i][offset_j+1]*rp2[i-1][j]);
            }
        }
        /************************Time step calculation - stop****************/
        /**********************Updating conservative and primitive variable - start ************/
        for (int j =1; j< jmax; j++)
        {
            for(int i =1; i< imax; i++)
            {
              ///conservative variables
                U[i][j][0] = U[i][j][0] - tvl[i][j]*Res[i][j][0];
                U[i][j][1] = U[i][j][1] - tvl[i][j]*Res[i][j][1];
                U[i][j][2] = U[i][j][2] - tvl[i][j]*Res[i][j][2];
                U[i][j][3] = U[i][j][3] - tvl[i][j]*Res[i][j][3];
                /// primitive
                q[i][j][3] = U[i][j][0];
                q[i][j][1] = U[i][j][1]/U[i][j][0];
                q[i][j][2] = U[i][j][2]/U[i][j][0];
                // q(:,:,1) = (rair-1)*(U(:,:,4)-0.5*(U(:,:,2).^2+U(:,:,3).^2)./U(:,:,1));
                q[i][j][0] = (r-1)*(U[i][j][3] - 0.5*(U[i][j][1]*U[i][j][1] + U[i][j][2]*U[i][j][2])/U[i][j][0]);
            }
        }
        /**********************Updating conservative and primitive variable - stop ************/

        /***********************Boundary conditions which includes MP - start **************/

        /// *************************left face - can be inflow or interface*********************
        if(left_rank == -2)
        {
            /// subsonic
            if(mInf < 1)
            {
                for(int j =0; j< jmax+1; j++)
                {
                    q[0][j][0] = q[1][j][0];// pressure is imposed from right
                    q[0][j][1] = uInf;
                    q[0][j][2] = vInf;
                    q[0][j][3] = dInf;

                    U[0][j][0] = q[0][j][3];
                    U[0][j][1] = q[0][j][3]*q[0][j][1];
                    U[0][j][2] = q[0][j][3]*q[0][j][2];
                    //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                    U[0][j][3] = q[0][j][0]/(r-1) + 0.5*(q[0][j][3]*(q[0][j][1]*q[0][j][1] + q[0][j][2]*q[0][j][2])) ;
                }
            }
            else
            {
            /// supersonic - implicitly satisfied as everything is initialized to infi values
               for(int j =0; j< jmax+1; j++)
                {
                    q[0][j][0] = pInf;// pressure is imposed from right
                    q[0][j][1] = uInf;
                    q[0][j][2] = vInf;
                    q[0][j][3] = dInf;

                    U[0][j][0] = q[0][j][3];
                    U[0][j][1] = q[0][j][3]*q[0][j][1];
                    U[0][j][2] = q[0][j][3]*q[0][j][2];
                    //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                    U[0][j][3] = q[0][j][0]/(r-1) + 0.5*(q[0][j][3]*(q[0][j][1]*q[0][j][1] + q[0][j][2]*q[0][j][2])) ;
                }
            }


        }
        else
        {
            /// message has to be passed and received
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int j =0; j< jmax+1; j++)
                {
                 Left_send[counter] = q[1][j][k];
                    counter++;
                }
            }
            /// send and then receive
            mpiSend(Left_send,left_len,mpiDType,left_rank,1);
            mpiRecv(Left_recv,left_len,mpiDType,left_rank,1);
            /// update primitive and conservative
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int j =0; j< jmax+1; j++)
                {
                    q[0][j][k] = Left_recv[counter];
                    counter++;
                }
            }
            /// updating conservative variables
            for(int j =0; j< jmax+1; j++)
            {
                U[0][j][0] = q[0][j][3];
                U[0][j][1] = q[0][j][3]*q[0][j][1];
                U[0][j][2] = q[0][j][3]*q[0][j][2];
                //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                U[0][j][3] = q[0][j][0]/(r-1) + 0.5*(q[0][j][3]*(q[0][j][1]*q[0][j][1] + q[0][j][2]*q[0][j][2])) ;
            }


        }
        /// **********************top face - can be boundary or interface*********************************************
        if(top_rank == -1)
        {
            // flow tangency
            for(int i =0; i<imax-1; i++)
            {
                offset_i = i +grid_bottom.x;
                offset_j = jmax-1;
                q[i+1][jmax][1] = q[i+1][jmax-1][1] - (2* (q[i+1][jmax-1][1]*Jn[offset_i][offset_j][0] + q[i+1][jmax-1][2]*Jn[offset_i][offset_j][1]))*Jn[offset_i][offset_j][0];
                q[i+1][jmax][2] = q[i+1][jmax-1][2] - (2* (q[i+1][jmax-1][1]*Jn[offset_i][offset_j][0] + q[i+1][jmax-1][2]*Jn[offset_i][offset_j][1]))*Jn[offset_i][offset_j][1];
            }

            for(int i =0; i < imax+1; i++)
            {
                q[i][jmax][0] =  q[i][jmax-1][0];// pressure
                q[i][jmax][3] =  q[i][jmax-1][3];// density
                /*
                    q(2:imax,jmax+1,2) = q(2:imax,jmax,2) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,1);
                    q(2:imax,jmax+1,3) = q(2:imax,jmax,3) - (2*(q(2:imax,jmax,2).*yarean(1:imax-1,jmax,1)+q(2:imax,jmax,3).*yarean(1:imax-1,jmax,2))).*yarean(1:imax-1,jmax,2);
                */
                //q[i+1][jmax][1] = q[i+1][jmax-1][1] - (2*(q[i+1][jmax-1][1])*Jn[][jmax-1])
                U[i][jmax][0] = q[i][jmax][3];
                U[i][jmax][1] =  q[i][jmax][3]*q[i][jmax][1];
                U[i][jmax][2] =  q[i][jmax][3]*q[i][jmax][2];
                U[i][jmax][3] = q[i][jmax][0]/(r-1) + 0.5*(q[i][jmax][3]*(q[i][jmax][1]*q[i][jmax][1] + q[i][jmax][2]*q[i][jmax][2]));
            }

        }
        else
        {
        // message passing
        counter = 0;
        for(int k =0; k<4; k++)
            {
                for(int i =0; i< imax+1; i++)
                {
                 Top_send[counter] = q[i][jmax-1][k];
                    counter++;
                }
            }
            /// send and then receive
            mpiSend(Top_send,top_len,mpiDType,top_rank,1);
            mpiRecv(Top_recv,top_len,mpiDType,top_rank,1);
            /// update primitive and conservative
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int i =0; i< imax+1; i++)
                {
                    q[i][jmax][k] = Top_recv[counter];
                    counter++;
                }
             }
            /// updating conservative variables
           for(int i =0; i < imax+1; i++)
            {
                U[i][jmax][0] = q[i][jmax][3];
                U[i][jmax][1] =  q[i][jmax][3]*q[i][jmax][1];
                U[i][jmax][2] =  q[i][jmax][3]*q[i][jmax][2];
                U[i][jmax][3] = q[i][jmax][0]/(r-1) + 0.5*(q[i][jmax][3]*(q[i][jmax][1]*q[i][jmax][1] + q[i][jmax][2]*q[i][jmax][2]));
            }


        }

        /// *******************************right face - can be outflow or interface****************************************
        if(right_rank == -3)
        {
            /// subsonic
            if(mInf < 1)
            {
                for(int j =0; j< jmax+1; j++)
                {
                    q[imax][j][0] = pInf;// pressure is imposed from right
                    q[imax][j][1] = q[imax-1][j][1];
                    q[imax][j][2] = q[imax-1][j][2];
                    q[imax][j][3] = q[imax-1][j][3];

                    U[imax][j][0] = q[imax][j][3];
                    U[imax][j][1] = q[imax][j][3]*q[imax][j][1];
                    U[imax][j][2] = q[imax][j][3]*q[imax][j][2];
                    //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                    U[imax][j][3] = q[imax][j][0]/(r-1) + 0.5*(q[imax][j][3]*(q[imax][j][1]*q[imax][j][1] + q[imax][j][2]*q[imax][j][2])) ;
                }
            }
            else
            {
                /// supersonic - implicitly satisfied as everything is initialized to infi values
                for(int j =0; j< jmax+1; j++)
                {
                    q[imax][j][0] = q[imax-1][j][0];
                    q[imax][j][1] = q[imax-1][j][1];
                    q[imax][j][2] = q[imax-1][j][2];
                    q[imax][j][3] = q[imax-1][j][3];

                    U[imax][j][0] = q[imax][j][3];
                    U[imax][j][1] = q[imax][j][3]*q[imax][j][1];
                    U[imax][j][2] = q[imax][j][3]*q[imax][j][2];
                    //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                    U[imax][j][3] = q[imax][j][0]/(r-1) + 0.5*(q[imax][j][3]*(q[imax][j][1]*q[imax][j][1] + q[imax][j][2]*q[imax][j][2])) ;
                }
            }
        }
        else
        {
          /// message has to be passed and received
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int j =0; j< jmax+1; j++)
                {
                 Right_send[counter] = q[imax-1][j][k];
                    counter++;
                }
            }
            /// receive and then send
            mpiRecv(Right_recv,left_len,mpiDType,right_rank,1);
            mpiSend(Right_send,left_len,mpiDType,right_rank,1);
            /// update primitive and conservative
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int j =0; j< jmax+1; j++)
                {
                    q[imax][j][k] = Right_recv[counter];
                    counter++;
                }
             }
            /// updating conservative variables
            for(int j =0; j< jmax+1; j++)
            {
                U[imax][j][0] = q[imax][j][3];
                U[imax][j][1] = q[imax][j][3]*q[imax][j][1];
                U[imax][j][2] = q[imax][j][3]*q[imax][j][2];
                //U(1,:,4) = (q(1,:,1)/(rair-1))+0.5*(q(1,:,4).*(q(1,:,2).^2+q(1,:,3).^2));
                U[imax][j][3] = q[imax][j][0]/(r-1) + 0.5*(q[imax][j][3]*(q[imax][j][1]*q[imax][j][1] + q[imax][j][2]*q[imax][j][2])) ;
            }

        }
        /// ***************************bottom face - can be boundary or interface*************************************
        if(bottom_rank == -1)
        {
            // flow tangency
            for(int i =0; i<imax-1; i++)
            {
                offset_i = i +grid_bottom.x;
                offset_j = 0;
                q[i+1][0][1] = q[i+1][1][1] - (2* (q[i+1][1][1]*Jn[offset_i][offset_j][0] + q[i+1][1][2]*Jn[offset_i][offset_j][1]))*Jn[offset_i][offset_j][0];
                q[i+1][0][2] = q[i+1][1][2] - (2* (q[i+1][1][1]*Jn[offset_i][offset_j][0] + q[i+1][1][2]*Jn[offset_i][offset_j][1]))*Jn[offset_i][offset_j][1];
            }

            for(int i =0; i < imax+1; i++)
            {
                q[i][0][0] =  q[i][1][0];// pressure
                q[i][0][3] =  q[i][1][3];// density
                U[i][0][0] = q[i][0][3];
                U[i][0][1] =  q[i][0][3]*q[i][0][1];
                U[i][0][2] =  q[i][0][3]*q[i][0][2];
                U[i][0][3] = q[i][0][0]/(r-1) + 0.5*(q[i][0][3]*(q[i][0][1]*q[i][0][1] + q[i][0][2]*q[i][0][2]));
            }
        }
        else
        {
            // message passing
            counter = 0;
            for(int k =0; k<4; k++)
            {
                for(int i =0; i< imax+1; i++)
                {
                    Bottom_send[counter] = q[i][1][k];
                    counter++;
                }
            }
            /// receive and send
            mpiRecv(Bottom_recv,top_len,mpiDType,bottom_rank,1);
            mpiSend(Bottom_send,top_len,mpiDType,bottom_rank,1);
            /// update primitive and conservative
            counter =0;
            for(int k =0; k<4; k++)
            {
                for(int i =0; i< imax+1; i++)
                {
                    q[i][0][k] = Bottom_recv[counter];
                    counter++;
                }
            }
            /// updating conservative variables
            for(int i =0; i < imax+1; i++)
            {
                U[i][0][0] = q[i][0][3];
                U[i][0][1] =  q[i][0][3]*q[i][0][1];
                U[i][0][2] =  q[i][0][3]*q[i][0][2];
                U[i][0][3] = q[i][0][0]/(r-1) + 0.5*(q[i][0][3]*(q[i][0][1]*q[i][0][1] + q[i][0][2]*q[i][0][2]));
            }
            //cout << "teja" << endl;

        }

        /***********************Boundary conditions which includes MP - stop **************/

    }
    /**************************************main iteration loop - ends ***************************/

    /**************************************storing data into text - starts***************************/
    dataType qNodes[imax][jmax][4];
    for(int j =0; j < jmax; j++)
    {
        for(int i =0; i< imax; i++)
        {
            for(int k =0; k< 4; k++)
                qNodes[i][j][k] = 0.25*(q[i][j][k] + q[i+1][j][k] + q[i][j+1][k] + q[i+1][j+1][k]);
        }
    }
    string str = "txt/pressure_" + patch::to_string(current_rank) + ".txt";
    char *c = &str[0u];
    std::fstream coord (c,fstream::out);
    if(coord.is_open())
    {
        for(int j = 0; j < jmax; j++)
        {
            for(int i = 0 ; i < imax; i++)
            {
                coord << qNodes[i][j][0] <<"  ";
            }
            coord << "\n";
        }
        coord.close();
    }
    //resnorm
    str = "txt/resnorm_" + patch::to_string(current_rank) + ".txt";
    c = &str[0u];
    std::fstream coord2 (c,fstream::out);
    if(coord2.is_open()){
    for(int i =0; i < in_flow.iterations;i++)
    {
    coord2 << resNorm[i] << " \n";

    }
    coord2.close();
    }
    /**************************************storing data into text - ends***************************/

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
    string strx = "txt/process_xcoord_" + patch::to_string(current_rank) + ".txt";
    string stry = "txt/process_ycoord_" + patch::to_string(current_rank) + ".txt";
    char *cx = &strx[0u];char *cy = &stry[0u];
    std::fstream coordx (cx,fstream::out);
    std::fstream coordy (cy,fstream::out);
    if(coordx.is_open() && coordy.is_open() )
    {
        for(int j = grid_bottom.y; j <= grid_top.y; j++)
        {
            for(int i = grid_bottom.x; i <= grid_top.x; i++)
            {
                coordx << x_cord[i][j] << "  ";
                coordy << y_cord[i][j] << "  ";
            }
            coordx << "\n";
            coordy << "\n";
        }
        coordx.close();
        coordy.close();
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
        grid_top.y = grid_bottom.y + row_points;
        grid_bottom.x = proc_column *column_points;
        grid_top.x = grid_bottom.x + column_points;
        if(proc_row == vert_proc-1)
        {
            grid_top.y = height - 1;
        }
        if(proc_column == horz_proc -1)
        {
            grid_top.x = width -1;
        }
 ///cout << "currnt process is " << current_rank << " with "<< grid_bottom.x << " " << grid_bottom.y << " "<< grid_top.x << " "<< grid_top.y << " " << endl;
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
        params.Set_params(0,0,2,2);
        int total_proc = params.horz_proc * params.vert_proc + 1; // one to handle boundary and other is parent process

        if (total_proc != size)
        {
            cout << "processes shoud be "<< total_proc << endl;
            return 0;
        }
        params.Load_grid("bumpgrid.txt");// laod grid
        params.Set_ranks(rank,size); // set ranks of surrounding processes
        params.Set_indexes(); // compute indexes of the process
        /*set inflow conditions*/
        params.in_flow.iterations = 20000;
        params.in_flow.cfl = 0.3;
        params.in_flow.press = 1013.26;
        params.in_flow.temp = 300;
        params.in_flow.R = 287;
        params.in_flow.r = 1.4;
        params.in_flow.mach = 1.4;
        params.in_flow.dens = params.in_flow.press/(params.in_flow.R*params.in_flow.temp);
        params.in_flow.u = params.in_flow.mach*sqrt(params.in_flow.r*params.in_flow.press/params.in_flow.dens);
        params.in_flow.v =0;
        params.in_flow.w =0;
        //cout << params.in_flow.dens << endl;
        params.Run_vanleer();
        params.Write_coord();
        //cout << MPI::TAG_UB <<endl;
    }
   // cout << "Process " << rank << " done and waiting!!" << endl;
    MPI::Finalize(); // end MPI
    return 0;
}
