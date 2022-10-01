/*
2D Lid- Driven Cavity using Finite Difference Method
*/

#include<stdio.h>
#include<math.h>
#define Nx 129    			  //Nodes in X Direction
#define Ny 129   			  //Nodes in Y Direction
#define err_max pow(10,-5)    // Max Error

//Starting Main Function
int main()
{
	int i,j,k,n=0;     // n is used as counter and i,j,k used as looping Variables
	
	double del_x,del_y,delsqr_x,delsqr_y,e,sum,c,beta,Re,U;  //e=error, beta=aspect ratio, Re=Reynolds No, U=Upper Lid Velocity
	
	double u[Nx+1][Ny+1],v[Nx+1][Ny+1],w_old[Nx+1][Ny+1],w_new[Nx+1][Ny+1],shi_old[Nx+1][Ny+1],shi_new[Nx+1][Ny+1];
	//u = velocity in X direction
	//v = velocity in Y direction
	//w_old & w_new = old and New values of Vorticity
	//shi_old & shi_new = old and new value of Stream Function
	
	del_x=1.0/(Nx-1);
	del_y=1.0/(Ny-1);
	delsqr_x=del_x*del_x;
	delsqr_y=del_y*del_y;
	beta=(del_x/del_y);
	Re=400;
	U=1;
	
	//Opening File for Stream Function & Velocity Component U,V
	FILE *u_comp,*v_comp,*stream_func,*w,*u_centre,*v_centre;
	stream_func=fopen("Stream Function.dat","w");
	u_comp=fopen("Velocity Components U.dat","w");
	v_comp=fopen("Velocity Components V.dat","w");
	w=fopen("Vorticity.dat","w");
	u_centre=fopen("Velocity Components U at Centre Line.dat","w");
	v_centre=fopen("Velocity Components V at Centre Line.dat","w");
	
	//Initailization Shi and W on all Entire Grid Nodes
	for(j=1; j<=Ny ; j++)
	{
		for(i=1 ; i<=Nx ; i++)
		{
			shi_old[i][j]=0;
			w_old[i][j]=0;	
		}		
	}
	
	for(i=2 ; i<Nx ; i++)
	{
		w_old[i][Ny]=-(2*U/del_y);
    }
	
	for(j=1; j<=Ny ; j++)
	{		
		for(i=1 ; i<=Nx ; i++)
		{
			shi_new[i][j]=shi_old[i][j];
			w_new[i][j]=w_old[i][j];
		}
	}
	
	//Initailization U and V on BC Grid Nodes at X=1 and X=L
	for(i=1 ; i<=Nx ; i++)
	{
		if(i==1)
		{
			for (j=1 ; j<=Ny ; j++)
			{
				u[i][j]=0;
				v[i][j]=0;
			}
		}
		else if(i==Nx)
		{
			for (j=1 ; j<Ny ; j++)
			{
				u[i][j]=0;
				v[i][j]=0;
			}
		}			
	}
	
	//Initailization U and V on BC Grid Nodes at Y=1 and Y=L
	for(j=1 ; j<=Ny ; j++)
	{
		if(j==1)
		{
			for (i=1 ; i<=Nx ; i++)
			{
				u[i][j]=0;
				v[i][j]=0;
			}
		}
		else if(j==Ny)
		{
			for (i=1 ; i<Nx ; i++)
			{
				u[i][j]=U;
				v[i][j]=0;
			}
		}			
	}
	
	//Start of Loop
	do
	{
		//Increase counter by 1
		n++;
		
		//Compute Shi
		for(k=1 ; k<=10 ; k++)
		{		
			for(j=2 ; j<Ny ; j++)
			{
				for(i=2 ; i<Nx ; i++)
				{
	      			shi_new[i][j]=((shi_new[i+1][j]+shi_new[i-1][j]+(beta*beta*(shi_new[i][j+1]+shi_new[i][j-1]))+(delsqr_x*w_new[i][j]))/(2*(1+beta*beta)));
				}
			}			
		}
		
		//Apply Boundary Conditions for vorticity
		for(j=2 ; j<Ny ; j++)
		{
			w_new[1][j]=-(2*shi_new[2][j])/delsqr_x;          			//left wall
			w_new[Nx][j]=-(2*shi_new[Nx-1][j])/delsqr_x;      			//Right wall
		}
		
		for(i=2 ; i<Nx ; i++)
		{
			w_new[i][1]=-(2*shi_new[i][2])/delsqr_y;                	//Bottom wall
			w_new[i][Ny]=-(2*shi_new[i][Ny-1]+2*del_y*U)/delsqr_y;    	//Top wall
		}
		
		//Compute Vorticity
		for(k=1 ; k<=2 ; k++)
		{	
			for(j=2 ; j<Ny ; j++)
			{
				for(i=2 ; i<Nx ; i++)
				{
          			w_new[i][j]=(w_new[i+1][j]+w_new[i-1][j]+beta*beta*(w_new[i][j+1]+w_new[i][j-1])-(0.25*beta*Re*(w_new[i+1][j]-w_new[i-1][j])*(shi_new[i][j+1]-shi_new[i][j-1]))+(0.25*beta*Re*(w_new[i][j+1]-w_new[i][j-1])*(shi_new[i+1][j]-shi_new[i-1][j])))/(2*(1+beta*beta));
				}
			}
		}
		
		//Error Calculation
		sum=0;
		c=0;
		
		for(j=1 ; j<=Ny ; j++)
		{
			for(i=1 ; i<=Nx ; i++)
			{
				sum=sum+fabs(w_new[i][j]-w_old[i][j]);
				c=c+fabs(w_new[i][j]);
			}
		}
		
		e=fabs(sum/c);
		printf("Error = %0.5lf\n",e);
		
		//Reassign Values
		for(j=1 ; j<=Ny ; j++)
		{
			for(i=1 ; i<=Nx ; i++)
			{
				shi_old[i][j]=shi_new[i][j];
				w_old[i][j]=w_new[i][j];
			}
		}	
				
	}
	while(e>err_max);
	
	//Print Stream Function Values in Output file
	for(j=1 ; j<=Ny ; j++)
	{
		for( i=1 ; i<=Nx ; i++)
		{
		 	fprintf(stream_func,"%d\t%d\t%lf\n",i,j,shi_new[i][j]);
		}
	}
	
	//Print Vorticity Values in Output file
	for(j=1 ; j<=Ny ; j++)
	{
		for( i=1 ; i<=Nx ; i++)
		{
		 	fprintf(w,"%d\t%d\t%lf\n",i,j,w_new[i][j]);
		}
	}

	//Computing U and V
	for(j=2 ; j<Ny ; j++)
	{
		for(i=2 ; i<Nx ; i++)
		{
		 	u[i][j]=(shi_new[i][j+1]-shi_new[i][j-1])/(2*del_y);
			v[i][j]=(-1)*(shi_new[i+1][j]-shi_new[i-1][j])/(2*del_x);
		}
	}

	//Printing U and V in Output file
	for(j=1 ; j<=Ny ; j++)
	{
		for(i=1 ; i<=Nx ; i++)
		{
			fprintf(u_comp,"%d\t%d\t%lf\n",i,j,u[i][j]);
			fprintf(v_comp,"%d\t%d\t%lf\n",i,j,v[i][j]);
        }
    }
    
    //Printing U in Output file at centre Line

	for(j=1 ; j<=Ny ; j++)
	{
		fprintf(u_centre,"%lf\t%d\n",u[65][j],j);
    }

    
    //Printing V in Output file at centre Line

	for(i=1 ; i<=Nx ; i++)
	{
		fprintf(v_centre,"%d\t%lf\n",i,v[i][65]);
    }


    //Counting No of Iterations
  	printf("\nTotal No of Iteration is : %d",n);
  	
  	//Closing Files
  	fclose(stream_func);
  	fclose(u_comp);
  	fclose(v_comp);
  	fclose(w);
	fclose(u_centre);
	fclose(v_centre);
	 	
	return 0;
	
}

