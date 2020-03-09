// criado em					: 2019/06/14
// ultima atualização	: 2019/07/18
// autor							: Mariana Casaroto <mariana.fcasaroto@gmail.com>
// notas							: Modelo de Lebwohl-Lasher com o método de Barker-Watts
// compilação					: gcc bw.c -lm -lgsl -lgslcblas
// execução						: ./a.out 'n'		-> em que n é um numero a ser utilizado como semente (opcional)


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>



const int     Nx =	  10;  // Tamanho da rede em x. 
const int     Ny =	  10;  // Tamanho da rede em y. 
const int     Nz =	  10;  // Tamanho da rede em z. 
const double  dT =  0.05;  // Variação de temperatura. 
const double  dTc=  0.01;  // Variação de temperatura perto da transição de fase. 
const double  Ti =	 0.1;  // Temperatura inicial. 
const double  Tf =   1.7;  // Temperatura final. 
const double  Tic=	1.05;  // Temperatura inicial perto da transição de fase. 
const double  Tfc=   1.2;  // Temperatura final perto da transição de fase. 
const int     MC =   1e4;  // Número de passos de monte-carlo. 
const int     NT =   1e4;  // Número de passos no transiente. 
const int     MCc=   4e4;  // Número de passos de monte-carlo perto da transição de fase. 
const int     NTc=   2e4;  // Número de passos no transiente perto da transição de fase. 
const double   e =   1.0;  // Constante eij. 
const double K_B =   1.0;  // Constante de Boltzmann. 


 
/*Definição do tipo de variável da rede.*/
typedef struct{
	double nx;
	double ny;
	double nz;
}n;
            
n *rede;
n n_r;

/* output da rede */
void	op_rede(double T,  n *rede){
	int i, j, k;
	FILE *file;
	char name[100];

	sprintf(name, "../data/rede/temperatura-%g.csv",	T*10);
 	file= fopen(name, "w");	
	fprintf(file, "i,j,k,nx,ny,nz \n");
	for(i= 0; i<Nx; i++){
		for(j= 0; j<Ny; j++){
			for(k= 0; k<Nz; k++){
						fprintf(file,"%d,%d,%d,%e,%e,%e\n",	i, j, k,rede[(i*Ny+j)*Nz+k].nx,	rede[(i*Ny+j)*Nz+k].ny,	rede[(i*Ny+j)*Nz+k].nz);
			}		
		}
	}
	fclose(file);
}

/* Definição do produto escalar */
double p_esc(n v1, n v2){
	return (v1.nx*v2.nx)+(v1.ny*v2.ny)+(v1.nz*v2.nz);
}

/* Definição da energia */
double E_lb(double vip, double vim, double vjp, double vjm, double vkp, double vkm){

	return	-e*((1.5*vip*vip-0.5)+(1.5*vim*vim-0.5)+
        		  (1.5*vjp*vjp-0.5)+(1.5*vjm*vjm-0.5)+
          	 	(1.5*vkp*vkp-0.5)+(1.5*vkm*vkm-0.5));
}

//-----------------------------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv){
		 int i, j, k, s, ii, jj, kk, ss;   													// Contadores.
		 int ip, im, jp, jm, kp, km;																// Vizinhos com condição periódica de contorno.
		 int eixo, acc, ref;																				// Eixo a ser sortedo, taxa de aceitação e de recusa.
		 int NTS, MCS;																							// Numero de passos no transiente e de monte carlo	
	double ang, theta, phi, theta_r, phi_r;												// Angulos das moleculas da rede e da rede alternativa.
	double dE, E_i, E_r; 																					// Variação de energia, Energia da rede, Energia da rede após o flip, 
	double E, E1, E2, E4;	    																		// Energia total, Energia Média, Média da energia ao quadrado e a quarta.
	double T, DT;         																				// Temperatura, variação da temperatura;
	double sig, binder;  																					// Variancia da energia e paramentro de binder;
	double S1, S2;																								// Parâmetro de Ordem e Parâmetro de Ordem ao quadrado;
	double nex, ney, nez;																					// Componentes do vetor após a rotação ao redor de um eixo específico		
	double dp_ip, dp_im, dp_jp, dp_jm, dp_kp, dp_km;							// Produtos escalares entre um ponto da rede os vizinhos + proximos.
	double dp_ip_r, dp_im_r, dp_jp_r, dp_jm_r, dp_kp_r, dp_km_r; 	// Produtos escalares entre um ponto da rede alternativa e os vizinhos + proximos.
	double Q_00, Q_01, Q_02, Q_10, Q_11, Q_12, Q_20, Q_21, Q_22; 	// Entradas da matriz do parâmetro de ordem do sistema


	/* Abertura dos arquivos */
	FILE *file1;	  
	char name1[60];
  sprintf(name1, "../data/output_%d.dat",  Nx);
  file1 = fopen (name1, "w");

  /*------------------------*/


	/* Gerador de numeros aleatórios */	
	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL); // Difinição da semente.
	gsl_rng *r= gsl_rng_alloc(gsl_rng_taus);
	rede= (n *	 ) calloc(Nx*Ny*Nz, sizeof(n	 ));		// Aloca memória da rede.
	/*------------------------------*/


	
	/* Inicio da rede */
	for(i= 0; i<Nx; i++){
		for(j= 0; j<Ny; j++){
			for(k= 0; k<Nz; k++){
				rede[(i*Ny+j)*Nz+k].nx= 0;
				rede[(i*Ny+j)*Nz+k].ny= 0;
				rede[(i*Ny+j)*Nz+k].nz= 1;
			} 
		}
	}	
	op_rede(0, rede);
	/*-----------------*/

	Q_00= Q_11=	Q_22=	Q_01=	Q_02=	Q_10=	Q_12=	Q_20=	Q_21= 0.0;		
	E1= E2= E4= E= E_i= E_r= 0.0;
	S1= S2= 0.0;
	acc= 0;
	ref= 0;
	binder= 0.0;


	/*-------Variação da Temperatura------*/
	for(T= Ti; T<Tf; T+= DT){
		if(T>=Tic && T<=Tfc){
			DT=dTc;
			NTS=NTc;
			MCS=MCc;
		}else{
			DT=dT;
			NTS=NT;
			MCS=MC;
		}

 	  /* Transiente */
		for(s= 0; s<NTS; s++){
			for(ss= 0; ss<Nx*Ny*Nz; ss++){
				// sorteio do ponto aleatório
    		i= gsl_rng_uniform(r)*Nx;
      	j= gsl_rng_uniform(r)*Ny;
      	k= gsl_rng_uniform(r)*Nz;

	     	/* Vizinhos na rede com condição periódica de contorno*/
      	ip= (i+1)%Nx;
      	im= (i-1+Nx)%Nx;
     		jp= (j+1)%Ny;
     		jm= (j-1+Ny)%Ny;
      	kp= (k+1)%Nz;
    	  km= (k-1+Nz)%Nz;
				

				/* Definição do vetor aleatório */
				eixo =(int) gsl_ran_flat (r, 1, 3.99999); //escolhe um eixo aleatorio x=1 y=2 z=3

				ang=(M_PI/2.0);
				if (((1.0*(Nx*Ny*Nz)-(1.0*ref))/(1.0*Nx*Ny*Nz))<0.5){ 
					ang=0.95*ang;  //desmarque se k=0 e k=lz-1 não se move
				}else{ 
					ang=1.05*ang;
				}
				
 	
				if (ang>M_PI/2){
			 		ang=M_PI/2;
				}
 

 				theta_r = gsl_ran_flat (r,-ang,ang);          //escolhe um angulo aleatorio entre -ang e ang
 				if (eixo==1){
					nex= rede[(i*Ny+j)*Nz+k].nx;
   				ney= rede[(i*Ny+j)*Nz+k].ny*cos(theta_r)+rede[(i*Ny+j)*Nz+k].nz*sin(theta_r);
   				nez=-rede[(i*Ny+j)*Nz+k].ny*sin(theta_r)+rede[(i*Ny+j)*Nz+k].nz*cos(theta_r);
				}else if (eixo==2){
 					nex= rede[(i*Ny+j)*Nz+k].nx*cos(theta_r)+rede[(i*Ny+j)*Nz+k].nz*sin(theta_r);
 					ney= rede[(i*Ny+j)*Nz+k].ny;
 					nez=-rede[(i*Ny+j)*Nz+k].nx*sin(theta_r)+rede[(i*Ny+j)*Nz+k].nz*cos(theta_r);	
   			}else{
     			nex= rede[(i*Ny+j)*Nz+k].nx*cos(theta_r)+rede[(i*Ny+j)*Nz+k].ny*sin(theta_r);
     			ney=-rede[(i*Ny+j)*Nz+k].nx*sin(theta_r)+rede[(i*Ny+j)*Nz+k].ny*cos(theta_r);
     			nez= rede[(i*Ny+j)*Nz+k].nz;
   			}

 				n_r.nx=nex;
 				n_r.ny=ney;
				n_r.nz=nez;

	
				dp_ip= p_esc(rede[(i*Ny+j)*Nz+k], rede[(ip*Ny+j)*Nz+k]);
				dp_im= p_esc(rede[(i*Ny+j)*Nz+k], rede[(im*Ny+j)*Nz+k]);
				dp_jp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jp)*Nz+k]);
				dp_jm= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jm)*Nz+k]);
				dp_kp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+kp]);
				dp_km= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+km]);

			
				dp_ip_r=  p_esc(n_r, rede[(ip*Ny+j)*Nz+k]);
				dp_im_r=  p_esc(n_r, rede[(im*Ny+j)*Nz+k]);
				dp_jp_r=  p_esc(n_r, rede[(i*Ny+jp)*Nz+k]);
				dp_jm_r=  p_esc(n_r, rede[(i*Ny+jm)*Nz+k]);
				dp_kp_r=  p_esc(n_r, rede[(i*Ny+j)*Nz+kp]);
				dp_km_r=  p_esc(n_r, rede[(i*Ny+j)*Nz+km]);
		
					
				E_i= E_lb(dp_ip, dp_im, dp_jp, dp_jm, dp_kp, dp_km);
				E_r= E_lb(dp_ip_r, dp_im_r, dp_jp_r, dp_jm_r, dp_kp_r, dp_km_r);
				
				dE=E_r-E_i;


				if(gsl_rng_uniform(r)< exp(-dE/(1.0*K_B*T))){
					rede[(i*Ny+j)*Nz+k].nx=n_r.nx;
					rede[(i*Ny+j)*Nz+k].ny=n_r.ny;
					rede[(i*Ny+j)*Nz+k].nz=n_r.nz;
					acc++; 
      	}else{
					ref++; 
				}
			}
		} /* Fim do transiente */

  	/* Calculo das energias iniciais */
		E=E_i=E_r= 0.0;
		binder= 0.0;
		Q_00= Q_11=	Q_22=	Q_01=	Q_02=	Q_10=	Q_12=	Q_20=	Q_21= 0.0;
		
    for(kk= 0; kk<Nz; kk++){
			for(ii= 0; ii< Nx; ii++){
				for(jj= 0; jj< Ny; jj++){

        	i= gsl_rng_uniform(r)*Nx;
        	j= gsl_rng_uniform(r)*Ny;
        	k= gsl_rng_uniform(r)*Nz;

        	// Vizinhos na rede com condição periódica de contorno
        	ip= (i+1)%Nx;
        	im= (i-1+Nx)%Nx;
        	jp= (j+1)%Ny;
        	jm= (j-1+Ny)%Ny;
        	kp= (k+1)%Nz;
        	km= (k-1+Nz)%Nz;

					dp_ip= p_esc(rede[(i*Ny+j)*Nz+k], rede[(ip*Ny+j)*Nz+k]);
					dp_im= p_esc(rede[(i*Ny+j)*Nz+k], rede[(im*Ny+j)*Nz+k]);
					dp_jp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jp)*Nz+k]);
					dp_jm= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jm)*Nz+k]);
					dp_kp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+kp]);
					dp_km= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+km]);


					E+= E_lb(dp_ip, dp_im, dp_jp, dp_jm, dp_kp, dp_km);

    		}
			}
		}


     /* Inicio do algorítmo de Monte Carlo*/
    for(s= 0; s< MCS; s++){
			for(ss= 0; ss<Nx*Ny*Nz; ss++){

      	// sorteio do ponto aleatório
      	i=  gsl_rng_uniform(r)*Nx;
      	j=  gsl_rng_uniform(r)*Ny;
      	k=  gsl_rng_uniform(r)*Nz;
		

				// Vizinhos na rede com condição periódica de contorno
      	ip= (i+1)%Nx;
      	im= (i-1+Nx)%Nx;
      	jp= (j+1)%Ny;
      	jm= (j-1+Ny)%Ny;
     		kp= (k+1)%Nz;
      	km= (k-1+Nz)%Nz;
 
				eixo =(int) gsl_ran_flat (r, 1, 3.99999); //escolhe um eixo aleatorio x=1 y=2 z=3
				
				ang=(M_PI/2.0);	
				if (((1.0*(Nx*Ny*Nz)-(1.0*ref))/(1.0*Nx*Ny*Nz))<0.5){ 
					ang=0.95*ang;  //desmarque se k=0 e k=lz-1 não se move
				}else{ 
					ang=1.05*ang;
				}
				

  			if (ang>M_PI/2.0){   //limita o angulo em pi/2
				 ang=M_PI/2.0;
				}

   			theta_r = gsl_ran_flat (r,-ang, ang);          //escolhe um angulo aleatorio entre -ang e ang
   			if (eixo==1){
					nex= rede[(i*Ny+j)*Nz+k].nx;
     			ney= rede[(i*Ny+j)*Nz+k].ny*cos(theta_r)+rede[(i*Ny+j)*Nz+k].nz*sin(theta_r);
     			nez=-rede[(i*Ny+j)*Nz+k].ny*sin(theta_r)+rede[(i*Ny+j)*Nz+k].nz*cos(theta_r);
				}else if (eixo==2){
     				nex= rede[(i*Ny+j)*Nz+k].nx*cos(theta_r)+rede[(i*Ny+j)*Nz+k].nz*sin(theta_r);
     				ney= rede[(i*Ny+j)*Nz+k].ny;
     				nez=-rede[(i*Ny+j)*Nz+k].nx*sin(theta_r)+rede[(i*Ny+j)*Nz+k].nz*cos(theta_r);	
   			}else{
     			nex= rede[(i*Ny+j)*Nz+k].nx*cos(theta_r)+rede[(i*Ny+j)*Nz+k].ny*sin(theta_r);
     			ney=-rede[(i*Ny+j)*Nz+k].nx*sin(theta_r)+rede[(i*Ny+j)*Nz+k].ny*cos(theta_r);
     			nez= rede[(i*Ny+j)*Nz+k].nz;
   			}

   			n_r.nx=nex;
   			n_r.ny=ney;
				n_r.nz=nez;

				
				// Produto escalar de um ponto aleátorio da rede com os vizinhos mais próximos
  	    dp_ip= p_esc(rede[(i*Ny+j)*Nz+k], rede[(ip*Ny+j)*Nz+k]);
    	  dp_im= p_esc(rede[(i*Ny+j)*Nz+k], rede[(im*Ny+j)*Nz+k]);
       	dp_jp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jp)*Nz+k]);
       	dp_jm= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+jm)*Nz+k]);
       	dp_kp= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+kp]);
       	dp_km= p_esc(rede[(i*Ny+j)*Nz+k], rede[(i*Ny+j)*Nz+km]);
 
 
        // Produto escalar do vetor aleatorio  com os vizinhos mais próximos de um ponto aleaório da rede
       	dp_ip_r=  p_esc(n_r, rede[(ip*Ny+j)*Nz+k]);
      	dp_im_r=  p_esc(n_r, rede[(im*Ny+j)*Nz+k]);
       	dp_jp_r=  p_esc(n_r, rede[(i*Ny+jp)*Nz+k]);
       	dp_jm_r=  p_esc(n_r, rede[(i*Ny+jm)*Nz+k]);
       	dp_kp_r=  p_esc(n_r, rede[(i*Ny+j)*Nz+kp]);
       	dp_km_r=  p_esc(n_r, rede[(i*Ny+j)*Nz+km]);
 
         
       	E_i= E_lb(dp_ip, dp_im, dp_jp, dp_jm, dp_kp, dp_km);
       	E_r= E_lb(dp_ip_r, dp_im_r, dp_jp_r, dp_jm_r, dp_kp_r, dp_km_r);
         
       	dE=E_r-E_i;
       	if(gsl_rng_uniform(r)< exp(-dE/(1.0*K_B*T))){
         rede[(i*Ny+j)*Nz+k].nx=n_r.nx;
         rede[(i*Ny+j)*Nz+k].ny=n_r.ny;
         rede[(i*Ny+j)*Nz+k].nz=n_r.nz;
				 acc++;  
      	}else{
       		ref++; 
       	}	

     		Q_00+= ((rede[(i*Ny+j)*Nz+k].nx)*(rede[(i*Ny+j)*Nz+k].nx)-(1.0/3.0));
     		Q_11+= ((rede[(i*Ny+j)*Nz+k].ny)*(rede[(i*Ny+j)*Nz+k].ny)-(1.0/3.0));
     		Q_22+= ((rede[(i*Ny+j)*Nz+k].nz)*(rede[(i*Ny+j)*Nz+k].nz)-(1.0/3.0));

     		Q_01+= (rede[(i*Ny+j)*Nz+k].nx)*(rede[(i*Ny+j)*Nz+k].ny);
     		Q_02+= (rede[(i*Ny+j)*Nz+k].nx)*(rede[(i*Ny+j)*Nz+k].nz);
     		Q_10+= (rede[(i*Ny+j)*Nz+k].ny)*(rede[(i*Ny+j)*Nz+k].nx);
     		Q_12+= (rede[(i*Ny+j)*Nz+k].ny)*(rede[(i*Ny+j)*Nz+k].nz);
     		Q_20+= (rede[(i*Ny+j)*Nz+k].nz)*(rede[(i*Ny+j)*Nz+k].nx);
     		Q_21+= (rede[(i*Ny+j)*Nz+k].nz)*(rede[(i*Ny+j)*Nz+k].ny);

			}/*Final de uma geração*/
			   
			/*Normalização das entradas do parâmetro de ordem por geração*/
			Q_00/= 1.0*(Nx*Ny*Nz);
    	Q_11/= 1.0*(Nx*Ny*Nz);
    	Q_22/= 1.0*(Nx*Ny*Nz);
    	Q_01/= 1.0*(Nx*Ny*Nz);
    	Q_02/= 1.0*(Nx*Ny*Nz);
    	Q_10/= 1.0*(Nx*Ny*Nz);
    	Q_12/= 1.0*(Nx*Ny*Nz);
    	Q_20/= 1.0*(Nx*Ny*Nz);
    	Q_21/= 1.0*(Nx*Ny*Nz);
		
			/* Calculo dos autovalores do Parâmetro de ordem tensorial */ 

    	double  data []= { Q_00, Q_01, Q_02,
      	                 Q_10, Q_11, Q_12,
        	               Q_20, Q_21, Q_22 };
	

    	gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
    	gsl_vector *eval = gsl_vector_alloc (3);                    // alloca memoria para os autovalores
    	gsl_matrix *evec = gsl_matrix_alloc (3, 3);                 // alloca memoria para os autovetores
    	gsl_eigen_symmv_workspace * r = gsl_eigen_symmv_alloc (3);  // alloca o espaço para computar os avt e avl  
    	gsl_eigen_symmv (&m.matrix, eval, evec, r);                 // função que computa os avt e avl 
    	gsl_eigen_symmv_free (r);                                   // libera memoria
    	gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_VAL_ASC);  // ordena os autovalores e autovetores

			/* Calculo do parâmetro de ordem escalar uniaxial S - maior autovalor*/
    	double eval_p = gsl_vector_get (eval, 2);
    	gsl_vector_view evec_p = gsl_matrix_column (evec, 2);
	
			S1+=eval_p;
			S2+=eval_p*eval_p;
	
	   	gsl_vector_free (eval);
			gsl_matrix_free (evec);
			/*-------------------------------------------------------*/

			/* Cálculo da energia após uma geração  */ 
			E=0;
			for(kk= 0; kk<Nz; kk++){
				for(ii= 0; ii< Nx; ii++){
					for(jj= 0; jj< Ny; jj++){

        		// Vizinhos na rede com condição periódica de contorno
        		ip= (ii+1)%Nx;
        		im= (ii-1+Nx)%Nx;
        		jp= (jj+1)%Ny;
        		jm= (jj-1+Ny)%Ny;
        		kp= (kk+1)%Nz;
        		km= (kk-1+Nz)%Nz;

						dp_ip= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(ip*Ny+jj)*Nz+kk]);
						dp_im= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(im*Ny+jj)*Nz+kk]);
						dp_jp= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(ii*Ny+jp)*Nz+kk]);
						dp_jm= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(ii*Ny+jm)*Nz+kk]);
						dp_kp= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(ii*Ny+jj)*Nz+kp]);
						dp_km= p_esc(rede[(ii*Ny+jj)*Nz+kk], rede[(ii*Ny+jj)*Nz+km]);


						E+=E_lb(dp_ip, dp_im, dp_jp, dp_jm, dp_kp, dp_km)/((double) Nx*Ny*Nz);
    			}
				}
			}	
	
			E1+=E;
			E2+=pow(E,2); 
			E4+=pow(E,4);
	
	
		}/* Final do Algorítmo de Metropolis*/
		
		
		op_rede(T, rede);		
 
		S1/=(1.0*MCS);
		S2/=(1.0*MCS);
		E1/=(1.0*MCS);
		E2/=(1.0*MCS); 	
		E4/=(1.0*MCS);
		sig=E2-(E1*E1);
		binder=1.0-(E4/(3.0*E2*E2));
		
    fprintf(file1,"%lf %lf %lf %lf %lf %lf %lf %lf %d %d\n", T, E1, E2, E4,(3.0/2.0)*S1, S2-S1*S1, sig, binder, acc, ref);
			
		S1=0.0;
			
	}/* Final da variação de Temperatura*/

	fclose(file1);
	free(rede);
	gsl_rng_free(r);

	return 0;
}




