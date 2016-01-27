#include <stdio.h>
#include <math.h>

// Defining speed of sound
double a(double T)
{
    double c;
    c = 1.4*287*T;
    return(c);
}

// Defining Mach
double Mach(double u,double v,double nx,double ny,double c)
{
    double M;
    M = (u*nx + v*ny)/c;
    return(M);
}

// Defining alpha +
double alp(double M)
{
    if (M > 0)
        return(1);
    else
        return(0);
}

// Defining alpha -
double aln(double M)
{
    if (M > 0)
        return(0);
    else
        return(1);
}

// Defining beta
double beta(double M)
{
    double c;
    c = abs(M);
    if(1-floor(c)>0)
        return(floor(c)-1);
    else
        return(0);
}

// Defining Cvlp
double Cvlp(double al,double be,double M)
{
    double cl;
    cl = al*(1+be)*M - (be/4)*pow(M+1,2) ;
    return(cl);
}

// Defining Cvlm
double Cvln(double al,double be,double M)
{
    double cl;
    cl = al*(1+be)*M + (be/4)*pow(1-M,2);
    return(cl);
}

// Defining Dbar +
double Dpos(double M,double be,double alp)
{
    double d;
    d = 0.25*(2-M)*pow(M+1,2);
    d = alp*(1+be) - be*d;
    return(d);
}

// Defining Dbar -
double Dneg(double M,double be,double aln)
{
    double d;
    d = 0.25*(2+M)*pow(M-1,2);
    d = aln*(1+be) - be*d;
    return(d);
}


int main()
{
    FILE *A = fopen("bumpgrid.txt","r");
    FILE *N = fopen("DATA_1_p1.txt","w"); // To store velocity
    FILE *O = fopen("Residue_1_p1.txt","w"); // To store residue
    FILE *B = fopen("Pressure_p1.txt","w");// To store the pressure
    int i,imx,jmx,j,k,n=0;
    double M_inf=0.4,u_inf,v_inf=0,Pin=1013.26,Tin=300,rho_inf=0.0117,cfl=0.5,res[25000]; // Given parameters
    double gamu,gamv,l,l1,lo,up,r,b,d,e;
    double Lx_1,Lx_2,Ly_1,Ly_2,Lamb_x1,Lamb_x2,Lamb_y1,Lamb_y2;
    fscanf(A,"%d %d",&imx,&jmx);
    u_inf = M_inf*sqrt(1.4*285*Tin);

    // Declaring u,v,fluxes,i-face normals,j-face normals,i-face areas,j-face areas,p,cell volumes..........
    double un[imx][jmx],vn[imx][jmx],pn[imx][jmx],Tn[imx][jmx],rhon[imx][jmx],M[imx][jmx]; // Node parameters
    double u[imx+1][jmx+1],v[imx+1][jmx+1],I[imx][jmx][2],J[imx][jmx][2],t[imx][jmx],p[imx+1][jmx+1],U[imx][jmx],Tem[imx+1][jmx+1],rho[imx+1][jmx+1],ho[imx+1][jmx+1];  // Cell centre parameters
    double c[imx][jmx],alphaxp[imx][jmx],alphaxn[imx][jmx],alphayp[imx][jmx],alphayn[imx][jmx],betx[imx][jmx],bety[imx][jmx]; // Speed of sound
    double Cvlxp[imx][jmx],Cvlxn[imx][jmx],Cvlyp[imx][jmx],Cvlyn[imx][jmx],Dxp[imx][jmx],Dxn[imx][jmx],Dyp[imx][jmx],Dyn[imx][jmx];
    double Q1[imx+1][jmx+1],Q2[imx+1][jmx+1],Q3[imx+1][jmx+1],Q4[imx+1][jmx+1];
    double F[imx][jmx][3],G[imx][jmx][3],Vol[imx][jmx],R[imx][jmx][3]; // i-face , j-face fluxes and Volume of the cells
    double Fp[imx][jmx][3],Fn[imx][jmx][3],Gp[imx][jmx][3],Gn[imx][jmx][3]; // Fluxes corresponding to positive U and negative U
    double Iar[imx][jmx],Jar[imx][jmx]; // i-face and j-face areas
    double X[imx][jmx],Y[imx][jmx]; // Position of the nodes

    // Initializing u,v and p on the entire domain (including ghost cells)
    for(j=0; j<jmx+1; j++)
    {
        for(i=0; i<imx+1; i++)
        {
            {
                p[i][j] = Pin;
                Q1[i][j] = rho_inf;
                Q2[i][j] = rho_inf*M_inf*sqrt(1.4*287*Tin);
                Q3[i][j] = rho_inf*v_inf;
                Q4[i][j] = rho_inf*Tin*(287/0.4) + 0.5*rho_inf*(u_inf*u_inf + v_inf*v_inf);
                rho[i][j] = Q1[i][j];
                u[i][j] = M_inf*sqrt(1.4*287*Tin);
                v[i][j] = 0;
                ho[i][j] = 3.5*(p[i][j]/rho[i][j]) + 0.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
            }

        }
    }

// Storing the coordinates
    for(j=0 ; j<jmx ; j++)
    {
        for(i=0 ; i<imx ; i++)
        {
            fscanf(A,"%lf %lf",&X[i][j],&Y[i][j]); // x and y coordinates
        }
    }

// Calculating i face areas
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx; i++)
        {
            Iar[i][j+1] = sqrt(pow((X[i][j+1] - X[i][j]),2)+ pow((Y[i][j+1] - Y[i][j]),2));
        }
    }

// Calculating j face areas
    for(j=0; j<jmx; j++)
    {
        for(i=0; i<imx-1; i++)
        {
            Jar[i+1][j] = sqrt(pow((X[i+1][j] - X[i][j]),2)+ pow((Y[i+1][j] - Y[i][j]),2));
        }
    }

// Calculating i face normals
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx; i++)
        {
            l = sqrt(pow((X[i][j+1] - X[i][j]),2)+ pow((Y[i][j+1] - Y[i][j]),2));
            I[i][j+1][0] = (Y[i][j+1] - Y[i][j])/l; // x-component of i-face normal
            I[i][j+1][1] = (-X[i][j+1] + X[i][j])/l; // y-component of i-face normal
        }
    }

// Calculating j face normals
    for(j=0; j<jmx; j++)
    {
        for(i=0; i<imx-1; i++)
        {
            l = sqrt(pow((X[i+1][j] - X[i][j]),2) + pow((Y[i+1][j] - Y[i][j]),2))	;
            J[i+1][j][0] = (-Y[i+1][j] + Y[i][j])/l; // x-component of j-face normal
            J[i+1][j][1] = (X[i+1][j] - X[i][j])/l; // y-component of j-face normal
        }
    }

    // Calculating cell volumes
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx-1; i++)
        {
            lo = 0.5*fabs((X[i+1][j] - X[i][j])*(Y[i+1][j+1] - Y[i+1][j]) - (X[i+1][j+1] - X[i+1][j])*(Y[i+1][j] - Y[i][j]));
            up = 0.5*fabs((X[i][j+1] - X[i][j])*(Y[i+1][j+1] - Y[i][j+1]) - (X[i+1][j+1] - X[i][j+1])*(Y[i][j+1] - Y[i][j]));
            Vol[i+1][j+1] = lo + up; // Volume of each cell
            printf("%f \n",Vol[i+1][j+1]);
        }
    }


    //Calculating c,M,alpha,beta,Cvl+,Cvl-......... for i faces
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx; i++)
        {
            c[i][j] = a(Tem[i][j]);
            M[i][j+1] = Mach(u[i][j],v[i][j],I[i][j+1][0],I[i][j+1][1],c[i][j]);
            alphaxp[i][j+1] = alp(M[i][j+1]);
            alphaxn[i][j+1] = aln(M[i][j+1]);
            betx[i][j+1] = beta(M[i][j+1]);
        }
    }
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx; i++)
        {
            Cvlxp[i][j+1] = Cvlp(alphaxp[i][j+1],betx[i][j+1],M[i][j+1]);
            Cvlxn[i][j+1] = Cvln(alphaxn[i+1][j+1],betx[i+1][j+1],M[i+1][j+1]);
            Dxp[i][j+1]=Dpos(M[i][j+1],betx[i][j+1],alphaxp[i][j+1]);
            Dxn[i][j+1]=Dneg(M[i+1][j+1],betx[i+1][j+1],alphaxn[i+1][j+1]);
            //printf("%lf %lf %d\n", Cvlxp[i][j+1], Cvlxn[i][j+1],i);
            printf("%d %d",i,j);
        }
    }

// Calculating c,M,alpha,beta,Cvl+,Cvl-......... for j faces
    for(i=0; i<imx-1; i++)
    {
        for(j=0; j<jmx; j++)
        {
            c[i][j] = a(Tem[i][j]);
            M[i+1][j] = Mach(u[i][j],v[i][j],J[i+1][j][0],J[i+1][j][1],c[i][j]);
            alphayp[i+1][j] = alp(M[i+1][j]);
            alphayn[i+1][j] = aln(M[i+1][j]);
            bety[i+1][j] = beta(M[i+1][j]);
        }
    }
    for(i=0; i<imx-1; i++)
    {
        for(j=0; j<jmx; j++)
        {
            Cvlyp[i+1][j] = Cvlp(alphayp[i+1][j],bety[i+1][j],M[i+1][j]);
            Cvlyn[i+1][j] = Cvln(alphayn[i+1][j+1],bety[i+1][j+1],M[i+1][j+1]);
            Dyp[i+1][j]=Dpos(M[i+1][j],bety[i+1][j],alphayp[i+1][j]);
            Dyn[i+1][j]=Dneg(M[i+1][j+1],bety[i+1][j+1],alphayn[i+1][j+1]);
        }
    }


// Flux reconsruction for i and j faces
    // Calculating F for i faces
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx; i++)
        {
            double qf[4]= {1,u[i][j+1],v[i][j+1],ho[i][j+1]} , qb[4]= {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]} , P[4] = {0,I[i][j+1][0],I[i][j+1][1],0};
            for (k=0; k<4; k++)
            {
                F[i][j+1][k] = rho[i][j+1]*c[i][j+1]*Cvlxp[i][j+1]*qf[k] + rho[i+1][j]*c[i+1][j]*Cvlxn[i+1][j+1]*qb[k] + (Dxp[i][j+1]*p[i][j+1] + Dxn[i+1][j+1]*p[i+1][j+1])*P[k] ;
            }
        }
    }

    // Calculating G for j faces
    for(i=0; i<imx-1; i++)
    {
        for(j=0; j<jmx; j++)
        {
            double qf[4]= {1,u[i+1][j],v[i+1][j],ho[i+1][j]} , qb[4]= {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]} , P[4] = {0,J[i+1][j][0],J[i+1][j][1],0};
            for (k=0; k<4; k++)
            {
                G[i+1][j][k] = rho[i+1][j]*c[i+1][j]*Cvlyp[i+1][j+1]*qf[k] + rho[i+1][j+1]*c[i+1][j+1]*Cvlyn[i+1][j+1]*qb[k] + (Dyp[i+1][j]*p[i][j] + Dyn[i+1][j+1]*p[i+1][j])*P[k] ;
            }
        }
    }

// Estimating time step for each cell
    for(j=1; j<jmx; j++)
    {
        for(i=1; i<imx; i++)
        {
            Lx_1 =  u[i-1][j]*I[i-1][j][0] + v[i-1][j]*I[i-1][j][1]; // Normal velocity across (i-1/2) face
            Lx_2 = u[i][j]*I[i][j][0] + v[i][j]*I[i][j][1]; // Normal velocity across (i+1/2) face
            Ly_1 = u[i][j-1]*J[i][j-1][0] + v[i][j-1]*J[i][j-1][1]; // Normal velocity across (j-1/2) face
            Ly_2 = u[i][j]*J[i][j][0] + v[i][j]*J[i][j][1]; // Normal velocity across (j+1/2) face
            Lamb_x1 = (fabs(Lx_1) + c[i][j])/2; // Max. normal velocity across (i-1/2) face
            Lamb_x2 = (fabs(Lx_2) + c[i][j])/2; // Max. normal velocity across (i+1/2) face
            Lamb_y1 = (fabs(Ly_1) + c[i][j])/2; // Max. normal velocity across (j-1/2) face
            Lamb_y2 = (fabs(Ly_2) + c[i][j])/2; // Max. normal velocity across (j+1/2) face
            t[i][j] = (cfl*2*Vol[i][j])/(Lamb_x1*Iar[i-1][j] + Lamb_x2*Iar[i][j] + Lamb_y1*Jar[i][j-1] + Lamb_y2*Jar[i][j]);
        }
    }


// Initial R(i,j)
    for(j=0; j<jmx-1; j++)
    {
        for(i=0; i<imx-1; i++)
        {
            for(k=0; k<4; k++)
            {
                R[i+1][j+1][k] = F[i+1][j+1][k] - F[i][j+1][k] + G[i+1][j+1][k] - G[i+1][j][k]; // All the 4 components of R(i,j)
            }
        }
    }

// Main Loop
    do
    {
        n=n+1;
        res[n] = 0;
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx-1; i++)
            {
                Q1[i+1][j+1] = Q1[i+1][j+1] - (t[i+1][j+1]/Vol[i+1][j+1])*(R[i+1][j+1][0]);  // Updating pressure
                Q2[i+1][j+1] = Q2[i+1][j+1] - (t[i+1][j+1]/Vol[i+1][j+1])*(R[i+1][j+1][1]); // Updating u
                Q3[i+1][j+1] = Q3[i+1][j+1] - (t[i+1][j+1]/Vol[i+1][j+1])*(R[i+1][j+1][2]); // Updating v
                Q4[i+1][j+1] = Q4[i+1][j+1] - (t[i+1][j+1]/Vol[i+1][j+1])*(R[i+1][j+1][3]); // Updating Et
            }
        }

// Updating u,v,rho,T and P
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx-1; i++)
            {
                rho[i][j] = Q1[i+1][j+1];
                u[i][j] = Q2[i+1][j+1]/Q1[i+1][j+1];
                v[i][j] = Q3[i+1][j+1]/Q1[i+1][j+1];
                Tem[i][j] = ((Q4[i+1][j+1]/rho[i][j]) - 0.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]))/717.5;
                p[i][j] = rho[i][j]*287*Tem[i][j];
            }
        }

// Inflow conditions
        for(j=0; j<jmx+1; j++)
        {
            Tem[0][j] = Tem[1][j];
            p[0][j] = rho[0][j]*287*Tem[0][j];
        }

// Outflow Conditions
        for(j=0; j<jmx+1; j++)
        {
            u[imx][j] = u[imx-1][j];
            v[imx][j] = v[imx-1][j];
            rho[imx][j] = rho[imx-1][j];
            p[imx][j] =  rho[imx][j]*287*Tem[imx][j];
        }

// Wall conditions
        for(i=0; i<imx+1; i++)
        {
            p[i][jmx] = p[i][jmx-1];
            p[i][0] = p[i][1];
            rho[i][jmx] = p[i][jmx]/(287*Tem[i][jmx]);
            rho[i][0] = p[i][0]/(287*Tem[i][0]);
        }

        // Calculating c,M,alpha,beta,Cvl+,Cvl-......... for i faces
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx; i++)
            {
                c[i][j] = a(Tem[i][j]);
                M[i][j+1] = Mach(u[i][j],v[i][j],I[i][j+1][0],I[i][j+1][1],c[i][j]);
                alphaxp[i][j+1] = alp(M[i][j+1]);
                alphaxn[i][j+1] = aln(M[i][j+1]);
                betx[i][j+1] = beta(M[i][j+1]);
            }
        }
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx; i++)
            {
                Cvlxp[i][j+1] = Cvlp(alphaxp[i][j+1],betx[i][j+1],M[i][j+1]);
                Cvlxn[i][j+1] = Cvln(alphaxn[i+1][j+1],betx[i+1][j+1],M[i+1][j+1]);
                Dxp[i][j+1]=Dpos(M[i][j+1],betx[i][j+1],alphaxp[i][j+1]);
                Dxn[i][j+1]=Dneg(M[i+1][j+1],betx[i+1][j+1],alphaxn[i+1][j+1]);
            }
        }

// Calculating c,M,alpha,beta,Cvl+,Cvl-......... for j faces
        for(i=0; i<imx-1; i++)
        {
            for(j=0; j<jmx; j++)
            {
                c[i][j] = a(Tem[i][j]);
                M[i+1][j] = Mach(u[i][j],v[i][j],J[i+1][j][0],J[i+1][j][1],c[i][j]);
                alphayp[i+1][j] = alp(M[i+1][j]);
                alphayn[i+1][j] = aln(M[i+1][j]);
                bety[i+1][j] = beta(M[i+1][j]);
            }
        }
        for(i=0; i<imx-1; i++)
        {
            for(j=0; j<jmx; j++)
            {
                Cvlyp[i+1][j] = Cvlp(alphayp[i+1][j],bety[i+1][j],M[i+1][j]);
                Cvlyn[i+1][j] = Cvln(alphayn[i+1][j+1],bety[i+1][j+1],M[i+1][j+1]);
                Dyp[i+1][j]=Dpos(M[i+1][j],bety[i+1][j],alphayp[i+1][j]);
                Dyn[i+1][j]=Dneg(M[i+1][j+1],bety[i+1][j+1],alphayn[i+1][j+1]);
            }
        }


// Flux reconsruction for i and j faces
        // Calculating F for i faces
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx; i++)
            {
                double qf[4]= {1,u[i][j+1],v[i][j+1],ho[i][j+1]} , qb[4]= {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]} , P[4] = {0,I[i][j+1][0],I[i][j+1][1],0};
                for (k=0; k<4; k++)
                {
                    F[i][j+1][k] = rho[i][j+1]*c[i][j+1]*Cvlxp[i][j+1]*qf[k] + rho[i+1][j]*c[i+1][j]*Cvlxn[i+1][j+1]*qb[k] + (Dxp[i][j+1]*p[i][j+1] + Dxn[i+1][j+1]*p[i+1][j+1])*P[k] ;
                }
            }
        }

        // Calculating G for j faces
        for(i=0; i<imx-1; i++)
        {
            for(j=0; j<jmx; j++)
            {
                double qf[4]= {1,u[i+1][j],v[i+1][j],ho[i+1][j]} , qb[4]= {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]} , P[4] = {0,J[i+1][j][0],J[i+1][j][1],0};
                for (k=0; k<4; k++)
                {
                    G[i+1][j][k] = rho[i+1][j]*c[i+1][j]*Cvlyp[i+1][j+1]*qf[k] + rho[i+1][j+1]*c[i+1][j+1]*Cvlyn[i+1][j+1]*qb[k] + (Dyp[i+1][j]*p[i][j] + Dyn[i+1][j+1]*p[i+1][j])*P[k] ;
                }
            }
        }

// Estimating time step for each cell
        for(j=1; j<jmx; j++)
        {
            for(i=1; i<imx; i++)
            {
                Lx_1 =  u[i-1][j]*I[i-1][j][0] + v[i-1][j]*I[i-1][j][1]; // Normal velocity across (i-1/2) face
                Lx_2 = u[i][j]*I[i][j][0] + v[i][j]*I[i][j][1]; // Normal velocity across (i+1/2) face
                Ly_1 = u[i][j-1]*J[i][j-1][0] + v[i][j-1]*J[i][j-1][1]; // Normal velocity across (j-1/2) face
                Ly_2 = u[i][j]*J[i][j][0] + v[i][j]*J[i][j][1]; // Normal velocity across (j+1/2) face
                Lamb_x1 = (fabs(Lx_1) + c[i][j])/2; // Max. normal velocity across (i-1/2) face
                Lamb_x2 = (fabs(Lx_2) + c[i][j])/2; // Max. normal velocity across (i+1/2) face
                Lamb_y1 = (fabs(Ly_1) + c[i][j])/2; // Max. normal velocity across (j-1/2) face
                Lamb_y2 = (fabs(Ly_2) + c[i][j])/2; // Max. normal velocity across (j+1/2) face
                t[i][j] = (cfl*2*Vol[i][j])/(Lamb_x1*Iar[i-1][j] + Lamb_x2*Iar[i][j] + Lamb_y1*Jar[i][j-1] + Lamb_y2*Jar[i][j]);
            }
        }


// Initial R(i,j)
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx-1; i++)
            {
                for(k=0; k<4; k++)
                {
                    R[i+1][j+1][k] = F[i+1][j+1][k] - F[i][j+1][k] + G[i+1][j+1][k] - G[i+1][j][k]; // All the 4 components of R(i,j)
                }
            }
        }

// Finding Residue
        for(j=0; j<jmx-1; j++)
        {
            for(i=0; i<imx-1; i++)
            {
                res[n] = res[n] +  pow((R[i+1][j+1][0]/(rho_inf*u_inf)),2) + pow((R[i+1][j+1][1]/(rho_inf*u_inf*u_inf)),2)+ pow((R[i+1][j+1][2]/(rho_inf*u_inf*u_inf)),2) + pow((R[i+1][j+1][3]/(rho_inf*u_inf*u_inf*u_inf)),2);
            }
        }
    }
//while(res[n]/res[1] > 0.00001 );
    while(n<10000);

    return(0);
}

