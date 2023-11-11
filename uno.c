//Programa que resuelve la ecuación de Schrödinger para el potencial V(x) de un pozo finito
//Encuentra el eigenvalor más pequeño usando el método de la potencia inversa 

#include <stdio.h>
#include<math.h>
#define NMAX 500


void LU ();

int main(){
double A[NMAX], C[NMAX], E[NMAX], u[NMAX];
double y[NMAX], q=100, q0, tol=1.0e-10, h, lt, e=4.7, V0=0.5 , V[NMAX];
double r=20.0; //buscamos integrar en [-r,r]
int n=400 , i , k=0;

h=2*r/(n+1);

FILE*arch=fopen("cinco.txt", "w");

lt=2*e;



for(i=0;i<n;i++){

    if(i>=(n/2)-r && i<=(n/2)+r){
        V[i]=0;
    }
    else{
        V[i]=V0;
    }

}

//  Aqui, graficamos el potencial para asegurarnos de que el condicional este definido correctamente


fprintf(arch, "%lf %lf \n", -r , V0 );

for (i = 0; i < n; i++) {
	fprintf(arch, "%lf %lf \n", -r+(i+1)*h , V[i] );
}

fprintf(arch, "%lf %lf \n", r , V0 );


/*
for(i=0;i<n;i++){
	
	A[i]=2.0*(1+(h*h)*V[i])/(h*h)-lt;
	E[i]=-1/(h*h);
	C[i]=-1/(h*h);
	u[i]=1/sqrt(n);
	
}

do {
	k++;
	LU(A,C,E,u,n,y);

	//normalizar a y
	q0=q;
	q=0;
	for (i = 0; i < n; i++) {
		q+=y[i]*y[i];
	}
	q=sqrt(q);

	for (i = 0; i < n; i++) {
		u[i]=y[i]/q;
	}


} while(fabs(q-q0)>tol);

printf("El eigenvalor más pequeño es: %lf \n", ((1.0/q)+lt)/2.0 );  //aquí cambiamos el signo dependiendo si nos acercamos por arriba o por abajo

fprintf(arch, "%lf %lf \n", -r , 0.0 );

for (i = 0; i < n; i++) {
	fprintf(arch, "%lf %lf \n", -r+(i+1)*h , u[i] );
}

fprintf(arch, "%lf %lf \n", r , 0.0 );

printf("El número de iteraciones es: %d\n", k);
*/
fclose(arch);

}




void LU (double A[NMAX], double C[NMAX], double E[NMAX], double b[NMAX], int n, double x[NMAX])
{
	double u[NMAX], w[NMAX], y[NMAX];
	int i;
	
	w[0]=A[0];
	
	for (i=0;i<n-1;i++){
		u[i]=E[i]/w[i];
		w[i+1]=A[i+1]-C[i+1]*u[i];
	}

	y[0]=b[0]/w[0];

	for (i=1;i<n;i++){
		y[i]=(b[i]-C[i]*y[i-1])/w[i];
	}

	x[n-1]=y[n-1];

	for (i=n-2;i>=0;i--){
		x[i]=y[i]-u[i]*x[i+1];
	}

}
