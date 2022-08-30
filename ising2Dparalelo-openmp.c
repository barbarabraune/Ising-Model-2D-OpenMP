#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

//definição de variaveis globais
#define L 500
#define IT 2000
#define EQ 100 // equilibra o sistema.

//inicializa o array 2D com uma configuração inicial de A
int **Inicializacao(int **A) 
{
  int i, j, x;
  srand(time(NULL) );
  for(i=0; i<L; i++)
  {
		for(j=0; j<L; j++)
    {
      x = -1+rand() % 3;

			if(x == 0)
      {
				A[i][j] = -1;
      }
      else 
				A[i][j] = 1;
    }
  }
	return A;
}

//calcula as condições de contorno
int **Contorno(int **A, int **B)
{
  int i, j;

  for(j=1; j<=L; j++)
  {
    B[0][j] = A[L-1][j-1];
    B[L+1][j] = A[0][j-1];
    B[j][0] = A[j-1][L-1];
    B[j][L+1] = A[L-1][0];
  }
  
  for(i=1; i<=L; i++)
  {
    for(j=1; j<=L; j++)
    {
      B[i][j] = A[i-1][j-1];
    }
  }

  return B;
}

//calcula a energia total da configuração de A
double Energia(int **B) 
{
  double en = 0, E = 0;
  int i, j;

  for(i=1; i<L+1; i++)
  {
		for(j=1; j<L+1; j++)
		{	
			en = B[i][j]*(B[i-1][j]+B[i][j+1]+B[i+1][j]+B[i][j-1]);
			E += -en;
		}
  }
	E = E/(L*L);		
  return E/2;
}

//calcula a magnetização por spins
double Magnetizacao(int **B)
{
  double total = 0;
  int i, j, N;

  N = L*L;

	for(i=1; i<L+1; i++)
  {
		for(j=1; j<L+1; j++)
    {
			total += B[i][j];
    }
  }
  return total/N;
}

void monteCarlo(int **B, double T, int chunk, int num_threads)
{
  int i, j, a, b;
  double dE, prob, random;
  omp_set_num_threads(num_threads);
    
  #pragma omp parallel private(i, j, a, b, dE, prob, random) shared(B, T) 
  {
    #pragma omp for schedule(dynamic, chunk) collapse(2)
  	for(i=1; i<L+1; i++)
  	{
  	  for(j=1; j<L+1; j++)
  	  {
  	    //flipar um spin aleatório
  	    a = 1+rand() % L;
  	    b = 1+rand() % L;
  	    //B[a][b] *= -1;

  	    //diferença de energia inicial e final
  	    dE = 2*B[a][b]*(B[a-1][b]+B[a][b+1]+B[a+1][b]+B[a][b-1]);

  	    //se delatE<=0, muda orientação
  	    if(dE<=0)
  	    {
  	      B[a][b] *= -1;
  	    }
  	    //senão, calcula prob=e^(-deltaE/T)
  	    else
  	    {
  	      prob = exp(-dE/T);
  	      //gera um número aleatório entre 0 e 1
  	      //random = -1+rand() % 3;
  	      random = drand48();

  	      //se o número aleatório for menor ou igual a probabilidade, muda a orientação do spin
  	      if(prob > random)
  	      {
  	        B[a][b]*=-1;
  	      }
  	    }
  	   }
  	  //fprintf (arq1, "%d %f %f\n", i, mag, Ef);
 	 }
   }

}

int main(int argc, char *argv[])
{
  FILE *arq1;
  int i, j, num_threads, chunk;
  int **A, **B, random, a, b, contagem=0;
  double prob, dE, Ei, Ef, mag, somaM=0, somaE=0, T;

  num_threads = atoi(argv[1]);
  chunk = IT/num_threads;

  arq1=fopen("grafico.dat", "wt");
  
  //aloca memória do array 2D
  A = malloc(L*sizeof(int*));
  for(i=0;i<L;i++) 
  {
    A[i] = (int*)malloc(L*sizeof(int));
  }

  //aloca memória do array 2D cópia
  B = malloc((L+2)*sizeof(int*));
  for(i=0;i<(L+2);i++) 
  {
    B[i] = (int*)malloc((L+2)*sizeof(int));
  }

  //chama para inicializar o array 2D
  A = Inicializacao(A);

  //aplica as condições de contorno
  B = Contorno(A, B);

  //calcula a magnetização inicial por spins
  mag = Magnetizacao(B);

  /*//imprime o array 2D para testar
  printf("\nConfiguração inicial dos spins:\n");
  for(int i=1; i<L+1; i++)
  {
    for(int j=1; j<L+1; j++)
    {
      printf("%d\t", B[i][j]);
    }
		printf ("\n");
  }*/

  //imprime a magnetização
  printf("\nMagnetização inicial: %.2lf\n", mag);

  //calcula a energia inicial 
  Ei = Energia(B);
  printf("\nEnergia inicial = %.2lf\n", Ei);
  
  /*printf("\nConfiguração depois da energia inicial:\n");
  for(int i=1; i<L+1; i++)
  {
    for(int j=1; j<L+1; j++)
    {
      printf("%d\t", B[i][j]);
    }
		printf ("\n");
  }*/
 
  srand(time(NULL) );
  for(T=0; T<4; T+=0.5)
  {
    A = Inicializacao(A);
    B = Contorno(A, B);
    
    for ( i = 0 ; i < EQ ; i++ )
    { // equilibra o sistema
      monteCarlo(B, T, chunk, num_threads);
    }

    mag = 0;
    Ef = 0;
    for ( i = 0 ; i < IT ; i++ )
    { 
      monteCarlo(B, T, chunk, num_threads);
      //soma a magnetização
      mag += Magnetizacao(B);
      //calcula a energia do novo estado
      Ef += Energia(B);
      //printf("\nMagnetização: %.2lf\n", mag);
    }
    
    fprintf (arq1, "%f %f %f\n", T, mag/IT, Ef/IT);
  }

  /*//imprime o array 2D
  printf("\nConfiguração final dos spins:\n");
  for(int i=1; i<L+1; i++)
  {
    for(int j=1; j<L+1; j++)
    {

      printf("%d\t", B[i][j]);
    }
		printf ("\n");
  }*/
  
  //imprime a magnetização final por spins
  printf("\nMagnetização final: %.2lf\n", mag);

  //energia final
  printf("\nEnergia final = %.2lf\n", Ef);

	return 0;
}
