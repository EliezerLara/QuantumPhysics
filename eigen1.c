//Programa que encuentra eigenvalores y eigenvectores mediante el método de potencias inversas

#include<stdio.h>
#include<math.h>

#define NMAX 100

void LU ();

//Función principal
int main () {

    double A[NMAX]={2,0} , C[NMAX]={0,1} , E[NMAX]={1,0} ;
    double u[NMAX]={1/sqrt(2),1/sqrt(2)}, y[NMAX], q=100 , q0 , eps=1.0e-10;
    int n=2 , i , k=0;

    do{
        LU(A,C,E,u,n,y);        //llamamos a la función LU
        
        q0=q;

            //Normalizamos a y
            q=0;

            for( i=0 ; i<n ; i++ ){

                q += y[i]*y[i];

            }

            q=sqrt(q);

            //aquí ya tenemos la normalizacion de y, que es q

            for( i=0 ; i<n ; i++ ){

                u[i]=y[i]/q;  //vector normalizado

            }

        
    k++;
    }
    while(fabs(q-q0)>eps);

    printf("El eigenvalor más pequeño es: %lf \n", 1.0/q);

    for( i=0; i<n ; i++ ){

        printf("%lf \n", u[i]);

    }

    printf("El número de iteraciones es: %d \n", k);

    return 0;
}



//Función  (subrutina)

void LU (double A[NMAX], double C[NMAX], double E[NMAX], double b[NMAX], int n , double x[NMAX]){
    
    double u[NMAX], w[NMAX], y[NMAX];

    int i;

    w[0]=A[0];

    for( i=0 ; i<n-1 ;i++){

        u[i]=E[i]/w[i];
        w[i+1]=A[i+1]-C[i+1]*u[i];

    }

    y[0]=b[0]/w[0];

    for( i=1; i<n ; i++ ){

        y[i]=(b[i]-C[i]*y[i-1])/w[i];
    }

    x[n-1]=y[n-1];

    for( i=n-2 ; i>=0 ; i-- ){

        x[i]=y[i]-u[i]*x[i+1];

    }


}