#include <stdio.h>
#include <math.h>

// Defining speed of sound
// double a(double p,double p1,double rho,double rho1)
 double a(double p,double rho)
  {
	double c;
	c = pow(1.4*p/rho,0.5);
	return(c);
  }


// Defining Mbar
 double Mbar(double M,double M1)
  {
	double Mb;
	Mb = sqrt(0.5*(pow(M,2)+pow(M1,2)));
	return(Mb);
  }

// Defining Mo
 double Mo(double M,double M1)
  {
	double M0;
	if(pow(M,2)>=pow(M1,2))
	  M0 = pow(M,2);
	   else
	    M0 = pow(M1,2);
	if(M0>=1)
	  M0 = 1;
	return(sqrt(M0));
  }

// Defining interface Mach for left running wave
 double Mach_L(double u,double v,double uf, double vf, double ub,double vb,double nx,double ny,double c)
  {
	double M,u_in,v_in,k = 1.0/3.0;
	u_in = u + (1-k)*(u-ub)/4.0 + (1+k)*(uf-u)/4.0;
	v_in = v + (1-k)*(v-vb)/4.0 + (1+k)*(vf-v)/4.0;
	M = (u_in*nx + v_in*ny)/c;
	return(M);
  }

// Define interface Mach for right running wave
double Mach_R(double u,double v,double uf, double vf, double ub,double vb,double nx,double ny,double c)
{
	double M,u_in,v_in,k = 1.0/3.0;
	u_in = u - (1-k)*(u-ub)/4.0 - (1+k)*(uf-u)/4.0;
	v_in = v - (1-k)*(v-vb)/4.0 - (1+k)*(vf-v)/4.0;
	M = (u_in*nx + v_in*ny)/c;
	return(M);
}
// ********************************************Defining Cx *********************************************
/*double C_x(double u,double v,double nx,double ny,double u1,double v1,double TO)
{
	double al,ar,M,M1,T_st,a_st;
	T_st = 2.0*TO/(2.4);
	a_st = sqrt(1.4*287*T_st);
	M = u*nx + v*ny;
	M1 = -(u1*nx + v1*ny);
	
	if(a_st >= M)
	   al = a_st;
    else
       al = pow(a_st,2)/M;
    
	if(a_st > M1)
	  ar = a_st;
	else
	  ar = pow(a_st,2)/(M1);
	
	if(al>=ar)
	  return(ar);
	else
	  return(al);   
}

//************************************************ Defining C- **************************************************
double C_y(double u,double v,double nx,double ny,double u1,double v1,double TO)
{
	double al,ar,M,M1,T_st,a_st;
	T_st = TO/(1.2);
	a_st = sqrt(1.4*287*T_st);
	M = u*nx + v*ny;
	M1 = -(u1*nx + v1*ny);
	
	if(a_st >= M)
	   al = a_st;
    else
       al = pow(a_st,2)/M;
    
	if(a_st > M1)
	  ar = a_st;
	else
	  ar = pow(a_st,2)/(M1);
	
	if(al>=ar)
	  return(ar);
	else
	  return(al);   
}

*/	
	

/*
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
*/

// ****************************************** Defining M+ *************************************************************
double Mp(double u,double v,double uf,double vf,double ub,double vb,double nx,double ny,double c)
{
	double M,Mplus;
	M = Mach_L(u,v,uf,vf,ub,vb,nx,ny,c);
	if(abs(M)>=1)
	   Mplus = 0.5*(M+abs(M));
	else
	   Mplus = 0.25*pow(M+1,2)*(0.5*pow(M-1,2) + 1);
	return(Mplus);
}

// ******************************************** Defining M- *************************************************************
double Mn(double u,double v,double uf,double vf,double ub,double vb,double nx,double ny,double c)
{
	double M,Mneg;
	M = Mach_R(u,v,uf,vf,ub,vb,nx,ny,c);
	if(abs(M)>=1)
	  Mneg = 0.5*(M-abs(M));
	else
	  Mneg = -0.25*pow(M-1,2)*(0.5*pow(M+1,2) + 1);
	return(Mneg); 
}

//******************************************* Defining M_half **********************************************************
double M_h(double u,double v,double nx,double ny,double uf,double vf,double uff,double vff,double ub,double vb,double c,double M_in,double rho,double rho1,double p,double p1)
{
	double M,M1,M_p,M_n,mbar,M_nt,fa,t,rh,M_half;
	M = Mach_L(u,v,uf,vf,ub,vb,nx,ny,c);
	M1 = Mach_R(uf,vf,uff,vff,u,v,nx,ny,c);
	M_p = Mp(u,v,uf,vf,ub,vb,nx,ny,c);
	M_n = Mn(uf,vf,uff,vff,u,v,nx,ny,c);
	mbar = Mbar(M,M1);
	M_nt = Mo(mbar,M_in);
	fa = M_nt*(2-M_nt);
	//fa=1;
	rh = 0.5*(rho+rho1);
	if((1-pow(mbar,2)) > 0)
	  t = (1-pow(mbar,2));
	else
	  t = 0;
	M_half = M_p + M_n - (0.25/fa)*(t)*(p1-p)/(rh*pow(c,2)) ;
	return(M_half);
}

// ****************************************** Defining P + ***************************************************************
 double Pp(double u,double v,double nx,double ny,double uf,double vf,double uff,double vff,double ub,double vb,double M_in,double c)
 {
 	double M,P_plus,fa,alpha,mbar,M1,M_nt;
 	M = Mach_L(u,v,uf,vf,ub,vb,nx,ny,c);
 	M1 = Mach_R(uf,vf,uff,vff,u,v,nx,ny,c);
 	mbar = Mbar(M,M1);
 	M_nt = Mo(mbar,M_in);
 	fa = M_nt*(2-M_nt);
 	//fa=1;
 	alpha = (3.0/16.0)*(-4.0+5.0*pow(fa,2));
 	if(abs(M)>=1)
 	  P_plus = (1.0/M)*(0.5*(M+abs(M)));
 	else
 	  P_plus = 0.25*pow(M+1,2)*((2-M) + 4.0*alpha*M*pow(M-1,2));
 	 // P_plus = 0.5*pow(M+1,3)*((3+M) + 4*alpha*M*pow(M-1,3)));
 	return(P_plus); 
 }
 
 // ***************************************** Defining P- ***************************************************************
 double Pn(double u,double v,double nx,double ny,double uf,double vf,double uff,double vff,double ub,double vb,double M_in,double c)
 {
 	double M,P_neg,fa,alpha,mbar,M1,M_nt;
 	M = Mach_L(u,v,uf,vf,ub,vb,nx,ny,c);
 	M1 = Mach_R(uf,vf,uff,vff,u,v,nx,ny,c);
 	mbar = Mbar(M,M1);
 	M_nt = Mo(mbar,M_in);
 	fa = M_nt*(2 - M_nt);
 	//fa=1;
 	alpha = (3.0/16.0)*(-4.0+5.0*pow(fa,2));	
 	if(abs(M1)>=1)
 	  P_neg = (1.0/M1)*(0.5*(M1-abs(M1)));
 	else
 	  P_neg = -0.25*pow(M1-1,2)*(-2-M1 + 4.0*alpha*M1*pow(M1+1,2));
 	return(P_neg);
 }


// ****************************************************  Defining Pu ***************************************************************
double Pu(double u, double v, double nx, double ny,double uf,double vf,double uff,double vff,double ub,double vb,double M_in,double c,double rho,double rho1,double p,double p1)
{
	double dpos,dneg,P,M,M1,mbar,fa,ul,ur,M_nt,P_p,P_n,ul_inf,ur_inf,vl_inf,vr_inf,k=1.0/3.0;
	P_p = Pp(u,v,nx,ny,uf,vf,uff,vff,ub,vb,M_in,c);
	P_n = Pn(u,v,nx,ny,uf,vf,uff,vff,ub,vb,M_in,c);
	ul_inf = u + (1-k)*(u-ub)/4.0 + (1+k)*(uf-u)/4.0;
	ur_inf = uf - (1-k)*(uf-u)/4.0 - (1-k)*(uff-uf)/4.0;
	vl_inf = v + (1-k)*(v-vb)/4.0 + (1+k)*(vf-v)/4.0;
	vr_inf = vf - (1-k)*(vf-v)/4.0 - (1-k)*(vff-vf)/4.0;
	ur = ur_inf*nx + vr_inf*ny;
	ul = ul_inf*nx + vl_inf*ny;
	M = Mach_L(u,v,uf,vf,ub,vb,nx,ny,c);
	M1 = Mach_R(uf,vf,uff,vff,u,v,nx,ny,c);
	mbar = Mbar(M,M1);
	M_nt = Mo(mbar,M_in);
	fa = M_nt*(2-M_nt);
	//fa=1;
	P =  P_p*p + P_n*p1 - 0.75*P_p*P_n*(rho + rho1)*fa*c*(ur-ul);
	return(P);
}


// ****************************************** Gf *****************************************************
double Gf(double q,double rho,double c,double cvl)
{
	double g;
	g=rho*c*cvl*q;
	return(g);
}

int main()
 {
	FILE *A = fopen("bumpgrid.txt", "r");
	FILE *Vel = fopen("Vel1.txt", "w"); // To store velocity
	FILE *Res = fopen("Residue1.txt", "w"); // To store residue
	FILE *Pre = fopen("Pressure1.txt", "w");// To store the pressure
	FILE *Mach = fopen("Mach1.txt","w"); //To store Mach
	int	i,j,k,n = 0;
	int imx,jmx;
	double M_in = 0.8, u_inf, v_inf = 0, Pin = 1013.26, Tin = 300, cfl=0.5, rho_inf, res[50100],T=0.000007; // Given parameters
	double To, l, l1, lo, up;
	double Lx_1, Lx_2, Ly_1, Ly_2, Lamb_x1, Lamb_x2, Lamb_y1, Lamb_y2;
	rho_inf = Pin/(287*Tin);
	fscanf(A, "%d %d", &imx, &jmx);
	u_inf = M_in*sqrt(1.4*287*Tin);
	printf("%lf %d %d \n",u_inf,imx,jmx);
    To = Tin*(1+0.2*pow(M_in,2));
    //printf("%lf \n",To);

//****************************************** Declaring u,v,fluxes,i-face normals,j-face normals,i-face areas,j-face areas,p,cell volumes.......... *******************************************************
	double un[imx][jmx],vn[imx][jmx],pn[imx][jmx],Mac[imx][jmx],rhon[imx][jmx]; // Node parameters
	double u[imx+1][jmx+1],v[imx+1][jmx+1],I[imx][jmx][2],J[imx][jmx][2],t[imx][jmx],p[imx+1][jmx+1],Px[imx][jmx],Py[imx][jmx],rho[imx+1][jmx+1],ho[imx+1][jmx+1];  // Cell centre parameters
	double c[imx+1][jmx+1],cx[imx][jmx],cy[imx][jmx]; // Speed of sound
//	double Cvlxp[imx][jmx],Cvlxn[imx][jmx],Cvlyp[imx][jmx],Cvlyn[imx][jmx],Dxp[imx][jmx],Dxn[imx][jmx],Dyp[imx][jmx],Dyn[imx][jmx];
    double Mx[imx][jmx],My[imx][jmx];
    double Q1[imx+1][jmx+1],Q2[imx+1][jmx+1],Q3[imx+1][jmx+1],Q4[imx+1][jmx+1];
	double F[imx][jmx][4],G[imx][jmx][4],Vol[imx][jmx],R[imx][jmx][4]; // i-face , j-face fluxes and Volume of the cells
	double Iar[imx][jmx],Jar[imx][jmx]; // i-face and j-face areas
	double X[imx][jmx],Y[imx][jmx]; // Position of the nodes


//*********************************************** Initializing u,v and p on the entire domain (including ghost cells) ***************************************************
	for(j=0; j<jmx+1; j++)
	{
	  for(i=0; i<imx+1; i++)
		{
			{
				p[i][j] = Pin;
				Q1[i][j] = rho_inf;
				Q2[i][j] = rho_inf*u_inf;
				Q3[i][j] = rho_inf*v_inf;
				Q4[i][j] = rho_inf*Tin*(287/0.4) + 0.5*rho_inf*(u_inf*u_inf + v_inf*v_inf);
				rho[i][j] = rho_inf;
				u[i][j] = u_inf;
				v[i][j] = v_inf;
				ho[i][j] = 3.5*(p[i][j]/rho[i][j]) + 0.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
				//printf("%d %d %lf \n",i,j,ho[i][j]);
			}

		}
	}

//********************************************************* Storing the coordinates ************************************************************************
	for(j=0; j<jmx; j++)
	{
	  for(i=0; i<imx; i++)
		{
			fscanf(A, "%lf %lf", &X[i][j], &Y[i][j]); // x and y coordinates
		}
	}

//********************************************************** Calculating i face areas ***********************************************************************
	for(j=0; j<jmx-1; j++)
	{
	  for(i=0; i<imx; i++)
		{
			Iar[i][j+1] = sqrt(pow((X[i][j+1] - X[i][j]), 2) + pow((Y[i][j+1] - Y[i][j]), 2));
		//	printf("%lf \n",Iar[i][j+1]);
		}
	}

//********************************************************** Calculating j face areas ************************************************************************
	for(j=0; j<jmx; j++)
	 {
		for(i=0; i<imx-1; i++)
		 {
			Jar[i+1][j] = sqrt(pow((X[i+1][j] - X[i][j]), 2) + pow((Y[i+1][j] - Y[i][j]), 2));
		 }
	 }

//******************************************************** Calculating i face normals *************************************************************************
	for(j=0; j<jmx-1; j++)
	 {
		for(i=0; i<imx; i++)
		 {
			l = sqrt(pow((X[i][j+1] - X[i][j]), 2) + pow((Y[i][j+1] - Y[i][j]), 2));
			I[i][j+1][0] = (Y[i][j+1] - Y[i][j])/l; // x-component of i-face normal
			I[i][j+1][1] = (-X[i][j+1] + X[i][j])/l; // y-component of i-face normal
			//printf("%lf \n",I[i][j + 1][1]);
		 }
	 }

//****************************************************** Calculating j face normals ***************************************************************************
	for(j=0; j<jmx; j++)
	 {
		for(i=0; i<imx-1; i++)
		 {
			l = sqrt(pow((X[i+1][j] - X[i][j]), 2) + pow((Y[i+1][j] - Y[i][j]), 2));
			J[i+1][j][0] = (-Y[i+1][j] + Y[i][j])/l; // x-component of j-face normal
			J[i+1][j][1] = (X[i+1][j] - X[i][j])/l; // y-component of j-face normal
		 }
	 }

//**************************************************** Calculating cell volumes ******************************************************************************
	for(j=0; j<jmx-1; j++)
	 {
		for(i=0; i<imx-1; i++)
		 {
			lo = 0.5*fabs((X[i+1][j] - X[i][j])*(Y[i+1][j+1] - Y[i+1][j]) - (X[i+1][j+1] - X[i+1][j])*(Y[i+1][j] - Y[i][j]));
			up = 0.5*fabs((X[i][j+1] - X[i][j])*(Y[i+1][j+1] - Y[i][j+1]) - (X[i+1][j+1] - X[i][j+1])*(Y[i][j+1] - Y[i][j]));
			Vol[i+1][j+1] = lo + up; // Volume of each cell
		 }
	 }

//************************************************** Finding c in each cell *********************************************************************************
	for (j = 0; j < jmx+1; j++)
	{
		for (i = 0; i < imx+1; i++)
		{
		 //	c[i][j] = a(p[i][j],p[i+1][j],rho[i][j],rho[i+1][j]);
		 c[i][j] = a(p[i][j],rho[i][j]);
		}
	}


//**************************************************** Calculation of Cvl+ , Cvl- , D+ ,D- for i faces ********************************************************
for(j=0; j<jmx-1; j++)
	 {
		for(i=0; i<imx; i++)
		 {
		//	k = 0.5*(c[i][j+1] + c[i+1][j+1]);
		//	Cvlxp[i][j + 1] = Cvlp(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,p[i][j+1],p[i+1][j+1],rho[i][j+1],rho[i+1][j+1],k);
		//	Cvlxn[i][j + 1] = Cvln(u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],k);
		//	Dxp[i][j + 1] = Dpos(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,k);
		//	Dxn[i][j + 1] = Dneg(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,k);
		    lo = p[i][j+1]-Pin;
		    up = p[i+1][j+1]-Pin;
	    	cx[i][j+1] = 0.5*(c[i][j+1] + c[i+1][j+1]);
		    if(i==0)
		     {
		       Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i][j+1],v[i][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i][j+1],v[i][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
	         }
			else if(i==imx-1)
			 {
	           Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i-1][j+1],v[i-1][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i-1][j+1],v[i-1][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
	         }
			else
			 {
	           Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i-1][j+1],v[i-1][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i-1][j+1],v[i-1][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
		     }
		//	 printf("%lf \n",cx[i][j+1]);
		 }
	 }

//************************************************** Calculation of Cvl+ , Cvl- , D+ ,D- for j faces ***********************************************************
	for(j=0; j<jmx; j++)
	{
		for(i=0; i<imx-1; i++)
		{
		//	k = 0.5*(c[i+1][j] + c[i+1][j+1]);
		//	Cvlyp[i + 1][j] = Cvlp(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,p[i+1][j],p[i+1][j+1],rho[i+1][j],rho[i+1][j+1],k);
		//	Cvlyn[i + 1][j] = Cvln(u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],k);
		//	Dyp[i + 1][j] = Dpos(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,k);
		//	Dyn[i + 1][j] = Dneg(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,k);
				lo = p[i+1][j]-Pin;
			 	up = p[i+1][j+1]-Pin;
			cy[i+1][j] = 0.5*(c[i+1][j] + c[i+1][j+1]);
			if (j==0)
			 {
			 	
			 	Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j],v[i+1][j],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j],v[i+1][j],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
			else if(j==jmx-1)
			 {
			  	Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i+1][j-1],v[i+1][j-1],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i+1][j-1],v[i+1][j-1],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
			else
			 {
				Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j-1],v[i+1][j-1],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j-1],v[i+1][j-1],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
		//	printf("%lf \n",cy[i+1][j]);
		}
	}


//*********************************************** Flux reconsruction for i and j faces ***********************************************************************

//*********************************************** Calculating G for j faces *********************************************************************************
	for(j=0; j<jmx ; j++)
	 {
		for(i=0; i<imx-1; i++)
		 {
			double m;
			double qf[4] = {1, u[i+1][j], v[i+1][j], ho[i+1][j] }, qb[4] = {1, u[i+1][j+1], v[i+1][j+1], ho[i+1][j+1]}, P[4] = { 0, J[i+1][j][0], J[i+1][j][1], 0 };
               // lo = 0.5*(c[i+1][j] + c[i+1][j+1]);
                
			    if (My[i+1][j] > 0)
			      m = cy[i+1][j]*rho[i+1][j]*My[i+1][j] ;
				else
				  m = cy[i+1][j]*rho[i+1][j+1]*My[i+1][j]; 

			for(k=0; k<4; k++)
		     {
				
                if(j==0 || j==jmx-1)
                 G[i+1][j][k] = Py[i+1][j]*P[k]*Jar[i+1][j];
				   else
                     if (m>0)
					   G[i+1][j][k] = m*qf[k]*Jar[i+1][j] + Py[i+1][j]*P[k]*Jar[i+1][j];
					    else
					     G[i+1][j][k] = m*qb[k]*Jar[i+1][j] + Py[i+1][j]*P[k]*Jar[i+1][j];
			   // printf("%d %d %d %lf \n",i,j,k,G[i+1][j][k]);
			 }
		 }
	}


//********************************************* Calculating F for i faces *************************************************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx; i++)
		   {
		   	double m;
		    double qf[4] = {1,u[i][j+1],v[i][j+1],ho[i][j+1]}, qb[4] = {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]}, P[4] = {0,I[i][j+1][0],I[i][j+1][1],0};
		     // lo = 0.5*(c[i][j+1] + c[i+1][j+1]);
			 
			   if(Mx[i][j+1] > 0)
			      m = cx[i][j+1]*rho[i][j+1]*Mx[i][j+1] ;
				else
				  m = cx[i][j+1]*rho[i+1][j+1]*Mx[i][j+1]; 
				  
		    for(k=0; k<4; k++)
			 {
			  if(m>0)
			   F[i][j+1][k] = m*qf[k]*Iar[i][j+1] + Px[i][j+1]*P[k]*Iar[i][j+1];
			    else
			      F[i][j+1][k] = m*qb[k]*Iar[i][j+1] + Px[i][j+1]*P[k]*Iar[i][j+1];
			// printf("%d %d %d %lf \n",i,j,k,F[i][j+1][k]);
			 }
		   }
		 }

//******************************************** Estimating time step for each cell *****************************************************************************
		for(j=1; j<jmx; j++)
		 {
		  for(i=1; i<imx; i++)
		   {
		   	//c[i][j] = 0.25*(cx[i-1][j] + cx[i][j] + cy[i][j-1] + cy[i][j]);
		    Lx_1 = u[i-1][j]*I[i-1][j][0] + v[i-1][j]*I[i-1][j][1]; // Normal velocity across (i-1/2) face
			Lx_2 = u[i][j]*I[i][j][0] + v[i][j]*I[i][j][1]; // Normal velocity across (i+1/2) face
			Ly_1 = u[i][j-1]*J[i][j-1][0] + v[i][j-1]*J[i][j - 1][1]; // Normal velocity across (j-1/2) face
			Ly_2 = u[i][j]*J[i][j][0] + v[i][j]*J[i][j][1]; // Normal velocity across (j+1/2) face
			Lamb_x1 = (fabs(Lx_1) + c[i][j]); // Max. normal velocity across (i-1/2) face
			Lamb_x2 = (fabs(Lx_2) + c[i][j]); // Max. normal velocity across (i+1/2) face
			Lamb_y1 = (fabs(Ly_1) + c[i][j]); // Max. normal velocity across (j-1/2) face
			Lamb_y2 = (fabs(Ly_2) + c[i][j]); // Max. normal velocity across (j+1/2) face
			t[i][j] = (cfl*2*Vol[i][j]) / (Lamb_x1*Iar[i-1][j] + Lamb_x2*Iar[i][j] + Lamb_y1*Jar[i][j-1] + Lamb_y2*Jar[i][j]);
		   }
		 }


//********************************************** Initial R(i,j) **************************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx-1; i++)
		   {
		    for(k=0; k<4; k++)
		     {
			  R[i+1][j+1][k] = F[i+1][j+1][k] - F[i][j+1][k] + G[i+1][j+1][k] - G[i+1][j][k]; // All the 4 components of R(i,j)
		     }
		  //  printf("%lf \n",R[i+1][j+1][1]);
		   }
		 }


//********************************************************* Main Loop ***********************************************************
   do
     {
		n = n + 1;
	    res[n] = 0;
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx-1; i++)
		   {
	     	Q1[i+1][j+1] = Q1[i+1][j+1] - (t[i+1][j+1]/ Vol[i+1][j+1])*(R[i+1][j+1][0]);  // Updating rho
			Q2[i+1][j+1] = Q2[i+1][j+1] - (t[i+1][j+1]/ Vol[i+1][j+1])*(R[i+1][j+1][1]); // Updating rho_u
			Q3[i+1][j+1] = Q3[i+1][j+1] - (t[i+1][j+1]/ Vol[i+1][j+1])*(R[i+1][j+1][2]); // Updating rho_v
			Q4[i+1][j+1] = Q4[i+1][j+1] - (t[i+1][j+1]/ Vol[i+1][j+1])*(R[i+1][j+1][3]); // Updating rho_Et
			
			// printf("%lf %lf %lf %lf \n",Q1[i + 1][j + 1],Q2[i + 1][j + 1],Q3[i + 1][j + 1],Q4[i + 1][j + 1]);
		   }
	     }

//************************************************* Updating u,v,rho,T and P ******************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx-1; i++)
		   {
		    rho[i+1][j+1] = Q1[i+1][j+1];
			u[i+1][j+1] = Q2[i+1][j+1]/Q1[i+1][j+1];
			v[i+1][j+1] = Q3[i+1][j+1]/Q1[i+1][j+1];
			p[i+1][j+1] = ((Q4[i+1][j+1]/rho[i+1][j+1]) - 0.5*(u[i+1][j+1]*u[i+1][j+1] + v[i+1][j+1]*v[i+1][j+1] ))*(0.4)*rho[i+1][j+1];
			//p[i+1][j+1] = 0.4*Q4[i+1][j+1];
			//if(n==30)
			//printf("%lf %lf %lf %lf \n",rho[i+1][j+1],u[i+1][j+1],v[i+1][j+1],p[i+1][j+1]);
     	   }
		 }

//****************************************************** Wall conditions ***********************************************************
		for (i=1; i<imx; i++)
		 {
		  p[i][jmx] = p[i][jmx-1];
		  p[i][0] = p[i][1];
		 }


// ****************************** Boundary conditions ************************************
      if(M_in<1.2)
        {
//****************************************************** Inflow conditions ****************************************************************
		  for (j=0; j<jmx+1; j++)
		   {
		     rho[0][j] = rho_inf;
		  	 p[0][j] = p[1][j];
		  //	 p[0][j] = Pin ;
		  	 u[0][j] = u_inf;
		  //	 u[0][j] = u[1][j];
		  	 v[0][j] = v_inf;
		   }

//****************************************************** Outflow Conditions ************************************************************
		  for (j=0; j<jmx+1; j++)
		   {
		     u[imx][j] = u[imx-1][j];
		   //  u[imx][j] = u_inf;
		  	// v[imx][j] = v_inf;
		  	   v[imx][j] = v[imx-1][j];
		  	// rho[imx][j] = rho_inf;
		  	 rho[imx][j] = rho[imx-1][j];
		  	 p[imx][j] = Pin;
		  	// p[imx][j] = p[imx-1][j];
	
		   }
		}
      else
        {
//****************************************************** Inflow conditions ****************************************************************
		  for (j=0; j<jmx+1; j++)
		   {
		     rho[0][j] = rho_inf;
		  	 p[0][j] = Pin;
		  	 u[0][j] = u_inf;
		  	 v[0][j] = v_inf;
		   }
        // printf("%lf \n",M_in);
//****************************************************** Outflow Conditions ************************************************************
		  for (j=0; j<jmx+1; j++)
		   {
		     u[imx][j] = u[imx-1][j];
		  	 v[imx][j] = v[imx-1][j];
		  	 rho[imx][j] = rho[imx-1][j];
		  	 p[imx][j] = p[imx-1][j];
		   }
		}

// ************************************ Boundary conditions for left running wave *************************************
/*		
		//****************************************************** Inflow conditions ****************************************************************
		for (j=0; j<jmx+1; j++)
		 {
		  rho[imx][j] = rho_inf;
		  p[imx][j] = p[imx-1][j];
		  u[imx][j] = u_inf;
		  v[imx][j] = v_inf;
		}

//****************************************************** Outflow Conditions ************************************************************
		for (j=0; j<jmx+1; j++)
		 {
		  u[0][j] = u[1][j];
		  v[imx][j] = v_inf;
		  rho[0][j] = rho[1][j];
		  p[0][j] = p[1][j];
		}
*/

	   /* for (j=0; j<jmx-1; j++)
		 {
		  for (i=0; i<imx-1; i++)
		   {
		   	//if(n==300)
		   //	printf("%lf %lf %lf %lf %d %d\n",rho[i + 1][j + 1],u[i + 1][j + 1],v[i + 1][j + 1],p[i + 1][j + 1],i,j);
		   }
	     }
*/	     

//********************************* Ghost cell velocities *************************************
            	u[0][0] = u[0][1];  // for inflow ghost cell boundary4
				v[0][0] = -v[0][1];
				u[0][jmx] = u[0][jmx-1];
				v[0][jmx] = -v[0][jmx-1];
				u[imx][0] = u[imx][1]; // for outflow ghost cell voundary
				v[imx][0] = -v[imx][1];
				u[imx][jmx] = u[imx][jmx-1];
				v[imx][jmx] = -v[imx][jmx-1];

				for(i=1; i<imx; i++)
				{
					l = u[i][1]*J[i][0][0] + v[i][1]*J[i][0][1];    // Flow tangency condition for bottom wall
					u[i][0] = u[i][1] - 2*l*J[i][0][0];
					v[i][0] = v[i][1] - 2*l*J[i][0][1];
					l1 = u[i][jmx-1]*J[i][jmx-1][0] + v[i][jmx-1]*J[i][jmx-1][1];  // Flow tangency condition for bottom wall
					u[i][jmx] = u[i][jmx-1] - 2*l1*J[i][jmx-1][0];
					v[i][jmx] = v[i][jmx-1] - 2*l1*J[i][jmx-1][1];
				}

//****************************************************** Finding ho **********************************************************************
		for (j=0; j<jmx+1; j++)
		 {
			for (i=0; i<imx+1; i++)
			 {
				ho[i][j] = 3.5*(p[i][j]/rho[i][j]) + 0.5*(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
			 }

		 }




//*************************************************** Finding c in each cell ***************************************************************
		for (j = 0; j<jmx+1; j++)
		 {
		  for (i = 0; i<imx+1; i++)
		   {
		    c[i][j] = pow(1.4*p[i][j]/rho[i][j],0.5);
		   // if(n==300)
		    // printf("%d %d %lf \n",i,j,c[i][j]);
		   }
		 }
    
    
//**************************************************** Calculation of Cvl+ , Cvl- , D+ ,D- for i faces ********************************************************
	for(j=0; j<jmx-1; j++)
	 {
		for(i=0; i<imx; i++)
		 {
		//	k = 0.5*(c[i][j+1] + c[i+1][j+1]);
		//	Cvlxp[i][j + 1] = Cvlp(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,p[i][j+1],p[i+1][j+1],rho[i][j+1],rho[i+1][j+1],k);
		//	Cvlxn[i][j + 1] = Cvln(u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],k);
		//	Dxp[i][j + 1] = Dpos(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,k);
		//	Dxn[i][j + 1] = Dneg(u[i][j+1],v[i][j+1],I[i][j + 1][0],I[i][j + 1][1],u[i+1][j+1],v[i+1][j+1],I[i][j + 1][0],I[i][j + 1][1],M_in,k);
		    lo = p[i][j+1]-Pin;
		    up = p[i+1][j+1]-Pin;
	    	cx[i][j+1] = 0.5*(c[i][j+1] + c[i+1][j+1]);
		    if(i==0)
		     {
		       Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i][j+1],v[i][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i][j+1],v[i][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
	         }
			else if(i==imx-1)
			 {
	           Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i-1][j+1],v[i-1][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i-1][j+1],v[i-1][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
	         }
			else
			 {
	           Px[i][j+1] = Pu(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i-1][j+1],v[i-1][j+1],M_in,cx[i][j+1],rho[i][j+1],rho[i+1][j+1],p[i][j+1],p[i+1][j+1]);
	           Mx[i][j+1] = M_h(u[i][j+1],v[i][j+1],I[i][j+1][0],I[i][j+1][1],u[i+1][j+1],v[i+1][j+1],u[i+2][j+1],v[i+2][j+1],u[i-1][j+1],v[i-1][j+1],cx[i][j+1],M_in,rho[i][j+1],rho[i+1][j+1],lo,up);
		     }
		//	 printf("%lf \n",cx[i][j+1]);
		 }
	 }

//************************************************** Calculation of Cvl+ , Cvl- , D+ ,D- for j faces ***********************************************************
	for(j=0; j<jmx; j++)
	{
		for(i=0; i<imx-1; i++)
		{
		//	k = 0.5*(c[i+1][j] + c[i+1][j+1]);
		//	Cvlyp[i + 1][j] = Cvlp(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,p[i+1][j],p[i+1][j+1],rho[i+1][j],rho[i+1][j+1],k);
		//	Cvlyn[i + 1][j] = Cvln(u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],k);
		//	Dyp[i + 1][j] = Dpos(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,k);
		//	Dyn[i + 1][j] = Dneg(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],J[i+1][j][0],J[i+1][j][1],M_in,k);
				lo = p[i+1][j]-Pin;
			 	up = p[i+1][j+1]-Pin;
			cy[i+1][j] = 0.5*(c[i+1][j] + c[i+1][j+1]);
			if (j==0)
			 {
			 	
			 	Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j],v[i+1][j],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j],v[i+1][j],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
			else if(j==jmx-1)
			 {
			  	Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i+1][j-1],v[i+1][j-1],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+1],v[i+1][j+1],u[i+1][j-1],v[i+1][j-1],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
			else
			 {
				Py[i+1][j] = Pu(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j-1],v[i+1][j-1],M_in,cy[i+1][j],rho[i+1][j],rho[i+1][j+1],p[i+1][j],p[i+1][j+1]);
			    My[i+1][j] = M_h(u[i+1][j],v[i+1][j],J[i+1][j][0],J[i+1][j][1],u[i+1][j+1],v[i+1][j+1],u[i+1][j+2],v[i+1][j+2],u[i+1][j-1],v[i+1][j-1],cy[i+1][j],M_in,rho[i+1][j],rho[i+1][j+1],lo,up);
		     }
		//	printf("%lf \n",cy[i+1][j]);
		}
	}


//*********************************************** Flux reconsruction for i and j faces ***********************************************************************

//*********************************************** Calculating G for j faces *********************************************************************************
	for (j=0; j<jmx ; j++)
	{
		for (i=0; i<imx-1; i++)
		{
			double m;
			double qf[4] = {1,u[i+1][j],v[i+1][j],ho[i+1][j]}, qb[4] = {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]}, P[4] = {0,J[i+1][j][0],J[i+1][j][1],0};
               // lo = 0.5*(c[i+1][j] + c[i+1][j+1]);
                
			    if (My[i+1][j] > 0)
			      m = cy[i+1][j]*rho[i+1][j]*My[i+1][j] ;
				else
				  m = cy[i+1][j]*rho[i+1][j+1]*My[i+1][j]; 

		if(j==0 || j==jmx-1)
		  {
		  	for (k=0; k<4; k++)
		  	  	G[i+1][j][k] = Py[i+1][j]*P[k]*Jar[i+1][j];
		  }
		  else
		    if(m>0)
		     {
		     	for (k=0; k<4; k++)
		     	 G[i+1][j][k] = m*qf[k]*Jar[i+1][j] + Py[i+1][j]*P[k]*Jar[i+1][j];
		     }
		      else
		       {
		       	 for (k=0; k<4; k++)
		       	  G[i+1][j][k] = m*qb[k]*Jar[i+1][j] + Py[i+1][j]*P[k]*Jar[i+1][j];
		       }
					   
			    //if(n==41)
			    //printf("%d %d %d %lf \n",i,j,k,G[i+1][j][k]);
		}
	}


//********************************************* Calculating F for i faces *************************************************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx; i++)
		   {
		   	double m;
		    double qf[4] = {1,u[i][j+1],v[i][j+1],ho[i][j+1]}, qb[4] = {1,u[i+1][j+1],v[i+1][j+1],ho[i+1][j+1]}, P[4] = {0,I[i][j+1][0],I[i][j+1][1],0};
		     // lo = 0.5*(c[i][j+1] + c[i+1][j+1]);
			 
			   if (Mx[i][j+1] > 0)
			      m = cx[i][j+1]*rho[i][j+1]*Mx[i][j+1] ;
				else
				  m = cx[i][j+1]*rho[i+1][j+1]*Mx[i][j+1]; 
				  
			if(m>0)
			 {
		   		for (k=0; k<4; k++)
		     	 F[i][j+1][k] = m*qf[k]*Iar[i][j+1] + Px[i][j+1]*P[k]*Iar[i][j+1];
		 	 }  
			else 
		 	 {
		 	 	 for(k=0;k<4;k++)
			      F[i][j+1][k] = m*qb[k]*Iar[i][j+1] + Px[i][j+1]*P[k]*Iar[i][j+1];
			 }
			/* for(k=0;k<4;k++)
			 {
			  if(n==8)
			  printf("%d %d %d %lf \n",i,j,k,F[i][j+1][k]);
		     }
		     */
		   }
		 }
		 
//******************************************** Estimating time step for each cell *****************************************************************************
		for(j=1; j<jmx; j++)
		 {
		  for(i=1; i<imx; i++)
	    	{
	    	  //c[i][j] = 0.25*(cx[i-1][j] + cx[i][j] + cy[i][j-1] + cy[i][j]);
			  Lx_1 = u[i-1][j]*I[i-1][j][0] + v[i-1][j]*I[i-1][j][1]; // Normal velocity across (i-1/2) face
			  Lx_2 = u[i][j]*I[i][j][0] + v[i][j]*I[i][j][1]; // Normal velocity across (i+1/2) face
			  Ly_1 = u[i][j-1]*J[i][j-1][0] + v[i][j-1]*J[i][j-1][1]; // Normal velocity across (j-1/2) face
			  Ly_2 = u[i][j]*J[i][j][0] + v[i][j]*J[i][j][1]; // Normal velocity across (j+1/2) face
			  Lamb_x1 = (fabs(Lx_1) + c[i][j]); // Max. normal velocity across (i-1/2) face
			  Lamb_x2 = (fabs(Lx_2) + c[i][j]); // Max. normal velocity across (i+1/2) face
		 	  Lamb_y1 = (fabs(Ly_1) + c[i][j]); // Max. normal velocity across (j-1/2) face
			  Lamb_y2 = (fabs(Ly_2) + c[i][j]); // Max. normal velocity across (j+1/2) face
			  t[i][j] = (cfl*2*Vol[i][j]) / (Lamb_x1*Iar[i-1][j] + Lamb_x2*Iar[i][j] + Lamb_y1*Jar[i][j-1] + Lamb_y2*Jar[i][j]);
			}
		 }

//********************************************** R(i,j) **************************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx-1; i++)
		  {
		   for(k=0; k<4; k++)
		    {
			 R[i+1][j+1][k] = F[i+1][j+1][k] - F[i][j+1][k] + G[i+1][j+1][k] - G[i+1][j][k]; // All the 4 components of R(i,j)
		    }
		  // if(n==4)
		 // printf("%lf %lf %lf %lf %d %d \n",R[i+1][j+1][0],R[i+1][j+1][1],R[i+1][j+1][2],R[i+1][j+1][3],i,j);
		  }
		}

//************************************************************ Residue Calculation ***********************************************************************
		for(j=0; j<jmx-1; j++)
		 {
		  for(i=0; i<imx-1; i++)
		   {
		    res[n] = res[n] + sqrt(pow((R[i+1][j+1][0] / (rho_inf*u_inf)), 2) + pow((R[i+1][j+1][1] / (rho_inf*u_inf*u_inf)), 2) + pow((R[i+1][j+1][2] / (rho_inf*u_inf*u_inf)), 2) + pow((R[i+1][j+1][3] / (rho_inf*u_inf*u_inf*u_inf)), 2));
		  }
		 }
		 printf("%d \n",n);
		 if(n%1000 == 0)
		   printf("%d %d \n",imx,jmx);
		}
 while(n < 33000);
 
  printf("%d %d \n",imx,jmx);

//***************************************** Nodal velocities estimation requires ghost cells velocities *********************************************************
//******************************* Boundary wall Ghost cell velocities are formulated below (including boundaries of inflow and outflow ghost cells) ***********************
	/*			u[0][0] = u[0][1];  // for inflow ghost cell boundary4
				v[0][0] = -v[0][1];
				u[0][jmx] = u[0][jmx-1];
				v[0][jmx] = -v[0][jmx-1];
				u[imx][0] = u[imx][1]; // for outflow ghost cell voundary
				v[imx][0] = -v[imx][1];
				u[imx][jmx] = u[imx][jmx-1];
				v[imx][jmx] = -v[imx][jmx-1];

				for(i=1; i<imx; i++)
				{
					l = u[i][1]*J[i][0][0] + v[i][1]*J[i][0][1];    // Flow tangency condition for bottom wall
					u[i][0] = u[i][1] - 2*l*J[i][0][0];
					v[i][0] = v[i][1] - 2*l*J[i][0][1];
					l1 = u[i][jmx-1]*J[i][jmx-1][0] + v[i][jmx-1]*J[i][jmx-1][1];  // Flow tangency condition for bottom wall
					u[i][jmx] = u[i][jmx-1] - 2*l1*J[i][jmx-1][0];
					v[i][jmx] = v[i][jmx-1] - 2*l1*J[i][jmx-1][1];
				}
    */
//***************************** Nodal velocities are determined by averaging the surrounding cells velocities ********************************************
				for(j=0; j<jmx; j++)
				{
					for(i=0; i<imx; i++)
					{
						un[i][j] = 0.25*(u[i][j] + u[i][j+1] + u[i+1][j] + u[i+1][j+1]);
						vn[i][j] = 0.25*(v[i][j] + v[i][j+1] + v[i+1][j] + v[i+1][j+1]);
						pn[i][j] = 0.25*(p[i][j] + p[i][j+1] + p[i+1][j] + p[i+1][j+1]);
						rhon[i][j] = 0.25*(rho[i][j] + rho[i][j+1] + rho[i+1][j] + rho[i+1][j+1]);
						//printf("%d %d %lf \n",i,j,p[i][j]);
						fprintf(Pre, "%d %d %lf\n", i, j, pn[i][j]);
					    fprintf(Vel, "%.9lf %.9lf %.9lf %.9lf\n", X[i][j], Y[i][j], un[i][j], vn[i][j]);
					}
				}

// ***************************************** Calculating and storing Mach number **********************************************
for(j=1;j<jmx;j++)
  {
 	for(i=1;i<imx;i++)
 	{
	   Mac[i][j] = pow(un[i][j]*un[i][j] + vn[i][j]*vn[i][j],0.5)/c[i][j];
 	   fprintf(Mach, "%.9lf %.9lf %.9lf %.9lf %.9lf \n", X[i][j], Y[i][j],c[i][j],rhon[i][j],Mac[i][j]);
    }
  }
//******************************************* Storing Residue ******************************************************************
				for(i=1; i<n; i++)
				{
						fprintf(Res, "%d %.9lf %.9lf %.9lf %.9lf %.9lf \n", i, log10(i), res[i]/res[1], log10(res[i]/res[1]),res[1],res[i]);
					//printf("%d %.9lf %.9lf %.9lf\n", i, log10(i), res[i] / res[1], log10(res[i] / res[1]));
				}

//********************************************** Closing files *****************************************************************
				fclose(Vel);
				fclose(A);
				fclose(Pre);
				fclose(Res);
  system("pause");
				return(0);
}





