/* 	Bayesian parameter estimation of FM modulated signals with arbitrary */
/* 	waveform modulation */
/*	
/*	13 January 2011	author: Gregory A. Kyriazis		 	   */

/*	Use Non Linear Least Squares to estimate the harmonic parameters
	of the modulating signal  */

/*	Find the maximum of the posterior density for the frequencies.  */
/*	Estimate the harmonic amplitudes. Evaluate the uncertainties    */
/*	associated with the estimates 								   */

/*	The software user supplies two things: the number of 	 	  */
/*	samples and an initial "starting guess" of the minimum (negative  */
/*	maximum) point, and values for the algorithm convergence 	   */
/*	parameters. Then the program searches for a local minimum,     */
/*	beginning from the starting guess, using the Direct Search 	  */
/*	algorithm of Hooke and Jeeves.  */

/*  This software was used in the report below for measuring
	FM modulated signals: */

/*	Kyriazis G. A., "Estimating parameters of complex modulating
	signals from prior information about their arbitrary waveform
	components," IEEE Trans. Instrum. Meas., Vol. 62, No. 6, 2013. */

/*  The author of this software is G.A. Kyriazis.		   */
/*  Permission to use, copy, modify, and distribute this software  */
/*  for any purpose without fee is hereby granted, provided that   */
/*  this entire notice is included in all copies of any software   */
/*  which is or includes a copy or modification of this software   */
/*  and in all copies of the supporting documentation for such	   */
/*  software.  					   */
					
/*	THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT    */
/*  ANY EXPRESS OR IMPLIED WARRANTY.  IN PARTICULAR, NEITHER THE   */
/*  AUTHOR NOR INMETRO MAKE ANY REPRESENTATION OR WARRANTY OF ANY	   */
/*  KIND CONCERNING THE MERCHANTABILITY OF THIS SOFTWARE OR ITS    */
/*  FITNESS FOR ANY PARTICULAR PURPOSE. 			   */
		

#include <analysis.h>
#include <cvirte.h>		/* Needed if linking in external compiler; harmless otherwise */
#include <userint.h>
#include <gpib.h>
#include <utility.h>   
#include <ansi_c.h>
#include <formatio.h>
#include <Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.h>

#define      VARS		250	/* max # of variables	     */
#define 	 N          31000
#define      NFUN       3
#define      TRUE       1
#define      FALSE      0

void f1(double);
void f2(double, double *);
void setGIJ_NLS (double);
void setGIJC (double, double *, double *);
void ampli_est (void);
void freq_ampli_est (void);
void ampli_NLS (void);
void result_freq_NLS (void);
void result_ampl_NLS (void);

double f(double *, int);
double best_nearby(double *, double *, double, int);
int hooke(int, double *, double *, double, double, int);

//void sampled_data (void); 
//void sampled_data1 (void); 
void GenerateWave (int);
void GenerateData (int);
void setGIJ (double *);
void prob (void);
void ortho (void);
void ampli (void);
void result_freq (void);
void result_ampl (void);
void ortho_err (void);

void Selecionar_menu (void);
void Configurar (void);
void Medir (void);
void Incerteza (void);
void Plotar (void);

int sis_painel;
int configuracao_painel;
int medicao_painel;
int incert_painel;
int graficos_painel;
int menu_barra;


/* global variables */

double	x_NLS;
double  omega2, omega_2, n_sigma;
int  	ifun, iters;
int 	shape;
double  data[N], data1[N], data_print[N], sinal[N], residue[N], *hi_NLS, *hiC;
double  *bi_NLS, *newbi, *new_bi, *end_bi;
double  **GIJ_NLS, **GIJC;
double  Freq_NLS, *Mod_NLS, *Fase_NLS, *Mod_corr_NLS, *Fase_corr_NLS;

double 	*GIJ_array, *TIJ_array, *M_array, *invM_array;
double 	*GIJC_array, *TIJC_array, *MC_array, *invMC_array; 
double 	*a_NLS, *b_NLS;

double  Acrms_NLS, THD, Mod_Fundamental_NLS;

int		funevals = 0;
double	x[VARS];
int		n;
double	delta[VARS], point[VARS];
double	prevbest;
int		nvars;
double	*startpt[VARS], *endpt[VARS];
int		nvars, itermax;
double	rho, epsilon;

double  dc, omega1, omega_1, phase2, phase_2, DeltaF, Delta_F;
double  zloge, h2bar, stloge, sigma, phat, st, nstloge, nst;
int     ino, nfun, jj;
int     novo;
double  Onda_Result[N], hi[NFUN];
double  D_vector[NFUN];
double  GIJ[N][NFUN], HIJ[N][NFUN], V_norm[NFUN][NFUN], invD_norm[NFUN][NFUN]; 
double  Modulo[(NFUN+1)/2], Fase[(NFUN+1)/2], Modulo_corr[(NFUN+1)/2], Fase_corr[(NFUN+1)/2];
double 	Dh2bar_djdk[3][3], V_eignorm[3][3];
double  D_eig[3], omega[3];


double  Acrms, Mod_Fundamental;

float 	Metcorr, Freqcorr, Range, Tsamp, Aper,
		Rsource, Cload, Df;
int 	i, j, k, n, p, s, t, q;
int 	Val[N];					 	

int 	novo, conta, count, objeto, item, ativo;
double  Freq_Gen, Freq_Gen_Mod, Fase_Gen, Fase_Gen_Mod, Ampl_Gen, Ampl_Gen_Mod, dc_Gen, dc_Gen_Mod;
int     Nharm_Gen_Mod;

double  Onda_Seno[N], Onda_Triang_Gen[N], Onda_Quad_Gen[N];    


main()

{

	sis_painel = LoadPanel (0, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", SYSPANEL);
	menu_barra = LoadMenuBarEx (sis_painel, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", MENU, __CVIUserHInst);
	configuracao_painel = LoadPanel (sis_painel, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", CONFIG);	
	medicao_painel = LoadPanel (sis_painel, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", MEDICAO);
	incert_painel = LoadPanel (sis_painel, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", INCERT);
	graficos_painel = LoadPanel (sis_painel, "Hooke_Jeeves_Student3_Wave_Fourier_Hz_FM_English.uir", GRAFICOS); 
	Selecionar_menu ();
	
}


void Selecionar_menu (void)

{
	int controle, objeto;

	novo = TRUE;
	DisplayPanel (sis_painel);

	while (TRUE) {
		if (novo)
		GetUserEvent (1, &objeto, &item);
			switch (item) {
				case MENU_CONFIGURACAO:
					Configurar ();
					break;
				case MENU_MEDICAO:
					Medir ();
					break;
				case MENU_INCERT:
					Incerteza ();
					break;
				case MENU_GRAFICOS:
					Plotar ();
					break;
				case MENU_SAIR:
					if (!novo) {
						novo = TRUE;
					}
					else {
						exit (1);
					}
					break;
			}
	}
}

//******************************CONFIGURAÇÃO DO MULTÍMETRO********************************************************************************


void Configurar (void) 

{

	int i;
	

	DisplayPanel (configuracao_painel);
	while (TRUE) {
 		GetUserEvent (1, &objeto, &item);
   		if (objeto == menu_barra) {
			switch (item) {
				case MENU_CONFIGURACAO:
					novo = FALSE;
					break;
				case MENU_MEDICAO:
					novo = FALSE;
					break;
				case MENU_INCERT:
					novo = FALSE;
					break;
				case MENU_GRAFICOS:
					novo = FALSE;
					break;
				case MENU_SAIR:
					novo = TRUE;
					HidePanel (configuracao_painel);
					HidePanel (medicao_painel);
					HidePanel (graficos_painel);
					HidePanel (incert_painel);
					break;
			}
		return;
		}
		if (objeto == configuracao_painel) {
			switch (item) {
				case CONFIG_CONFIGURAR:
					count += 1;
					if (count>1) {
						SetCtrlAttribute (medicao_painel, MEDICAO_MEDIR_FLICKER, ATTR_DIMMED, 1);
						SetCtrlAttribute (medicao_painel, MEDICAO_FREQ_NLS, ATTR_DIMMED, 1);
						SetCtrlAttribute (medicao_painel, MEDICAO_AMPL_NLS, ATTR_DIMMED, 1);
						SetCtrlAttribute (medicao_painel, MEDICAO_FREQ, ATTR_DIMMED, 1);
						SetCtrlAttribute (medicao_painel, MEDICAO_AMPL, ATTR_DIMMED, 1);
					}
					GetCtrlVal (configuracao_painel, CONFIG_INO, &ino);
					GetCtrlVal (configuracao_painel, CONFIG_FREQ_GEN, &Freq_Gen);
					GetCtrlVal (configuracao_painel, CONFIG_FREQ_GEN_MOD, &Freq_Gen_Mod);
					GetCtrlVal (configuracao_painel, CONFIG_FASE_GEN, &Fase_Gen);
					GetCtrlVal (configuracao_painel, CONFIG_FASE_GEN_MOD, &Fase_Gen_Mod);
					GetCtrlVal (configuracao_painel, CONFIG_AMPL_GEN, &Ampl_Gen);
					GetCtrlVal (configuracao_painel, CONFIG_AMPL_GEN_MOD, &Ampl_Gen_Mod);
					GetCtrlVal (configuracao_painel, CONFIG_DC_GEN, &dc_Gen);
					GetCtrlVal (configuracao_painel, CONFIG_DC_GEN_MOD, &dc_Gen_Mod);
					GetCtrlVal (configuracao_painel, CONFIG_NHARM_GEN_MOD, &Nharm_Gen_Mod);
	  				GetCtrlVal (configuracao_painel, CONFIG_SIGMA, &n_sigma);	  
					GetCtrlVal (configuracao_painel, CONFIG_DELTAF, &Delta_F);		// DESVIO DE FREQUÊNCIA GERADO
		            GetCtrlVal (configuracao_painel, CONFIG_WFM, &shape);
        		    GetCtrlVal (configuracao_painel, CONFIG_SAMPLE_SPACING, &Tsamp);
        		    switch (shape) {
                		case 0 :
                    		GenerateWave (1);
	                    	break;
                		case 1 :
                    		GenerateWave (2);
    	                	break;
                		case 2 :
                    		GenerateWave (3);
        	            	break;
	                }
					ifun = 2*Nharm_Gen_Mod+1;
					nfun = NFUN;
					SetCtrlAttribute (medicao_painel, MEDICAO_MEDIR, ATTR_DIMMED, 0);
					hi_NLS = (double *) malloc (ifun*sizeof(double));
					hiC = (double *) malloc ((ifun+1)*sizeof(double));
					Mod_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));
					Fase_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));
					Mod_corr_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));
					Fase_corr_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));
				    GIJ_NLS = (double **) malloc (ino*sizeof(double*));
	  				for (i=0; i<ino; i++) {
	  	 			   GIJ_NLS[i] = (double *) malloc (ifun*sizeof(double));
	  				}
				    GIJC = (double **) malloc (ino*sizeof(double*));
	  				for (i=0; i<ino; i++) {
	  	 			   GIJC[i] = (double *) malloc ((ifun+1)*sizeof(double));
	  				}
	  				M_array = (double *) malloc (ifun*ifun*sizeof(double));
					invM_array = (double *) malloc (ifun*ifun*sizeof(double));
	  				TIJ_array = (double *) malloc (ino*ifun*sizeof(double));
	  				GIJ_array = (double *) malloc (ino*ifun*sizeof(double));
			  	    MC_array = (double *) malloc ((ifun+1)*(ifun+1)*sizeof(double));
			  	    invMC_array = (double *) malloc ((ifun+1)*(ifun+1)*sizeof(double));
			  	    TIJC_array = (double *) malloc (ino*(ifun+1)*sizeof(double));
			  	    GIJC_array = (double *) malloc (ino*(ifun+1)*sizeof(double));
			  	    a_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));	  
	   				b_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));	  
					bi_NLS = (double *) malloc ((ifun+1)*sizeof(double));
					newbi = (double *) malloc (ifun*sizeof(double));
					new_bi = (double *) malloc ((ifun+1)*sizeof(double));
					end_bi = (double *) malloc ((ifun+1)*sizeof(double));
					break;
			}
		}
	}
	
}



//******************************GERAÇÃO DA FORMA DE ONDA***********************************************


void GenerateWave (int shape)
{

	int i, j, k;
	double Gaussrnd[N];
	double Gaussrnd_pha[N];
	double  Onda_Seno[N], Onda_Triang_Gen[N], Onda_Quad_Gen[N];    
    
    switch (shape) {
        case 1 :        /* Sinewave */
		    for (i=0; i<ino; i++) {
	  	 		data[i] = 0;
	  	 		Onda_Seno[i] = 0;
	 	  	}
	  		GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		for (i=0; i<ino; i++) {	  
	  	 		Onda_Seno[i] = Ampl_Gen_Mod*cos(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod);
			    Onda_Seno[i] += dc_Gen_Mod;
		  	}   
	  		for (i=0; i<ino; i++) {
		 		data[i] = Onda_Seno[i]+Gaussrnd[i];
	  		}
            break;
        case 2 :        /* Square wave */
		  	for (i=0; i<ino; i++) {
	  	 		data[i] = 0;
	  	 		Onda_Quad_Gen[i] = 0;
	  		}
	  		GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		for (i=0; i<ino; i++) {	  
  	 			for (j=0; j<Nharm_Gen_Mod/2; j++) {	  // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    		k = 2*j+1;
	  	    		Onda_Quad_Gen[i] += Ampl_Gen_Mod*sin(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod))/k;
	  	 		}
	     		Onda_Quad_Gen[i] *= 2/asin(1);
	     		Onda_Quad_Gen[i] += dc_Gen_Mod;
	  		}   
	  		for (i=0; i<ino; i++) {
		 		data[i] = Onda_Quad_Gen[i]+Gaussrnd[i];
	  		}
            break;
        case 3 :        /* Triangle wave */
		  	for (i=0; i<ino; i++) {
	  	 		data[i] = 0;
	  	 		Onda_Triang_Gen[i] = 0;
	  		}
	  		GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		for (i=0; i<ino; i++) {	  
	  	 		for (j=0; j<Nharm_Gen_Mod/2; j++) {	  // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    		k = 2*j+1;
	  	    		Onda_Triang_Gen[i] += Ampl_Gen_Mod*cos(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod))/pow(k, 2);
	  	    	}
	     		Onda_Triang_Gen[i] *= 2/pow(asin(1), 2);
	     		Onda_Triang_Gen[i] += dc_Gen_Mod;
	     	}   
	  		for (i=0; i<ino; i++) {
		 		data[i] = Onda_Triang_Gen[i]+Gaussrnd[i];
		 	}
            break;
    }        

}



//******************************MEDIÇÃO DA FREQUÊNCIA E DAS AMPLITUDES E FASES DOS HARMÔNICOS*************************************************************


void Medir (void)

{

	/*	VALOR EFICAZ ESPERADO DO SINAL A SER MEDIDO +/- 10%.
			A MENOR TENSÃO MENSURÁVEL É 1 mV, A MAIOR É 700 V.*/


	double	   startpt[VARS], endpt[VARS];
	double	   startpt_NLS;
    int		   itermax_NLS;
    double	   rho_NLS, epsilon_NLS;
    int		   vars, itermax;
    double	   rho, epsilon;
    double	   steplength_NLS;
    int		   j, k;
    double     Freq_start, Freq_start_carr, Fase_start, DeltaF_start;
    
   	char diretorio[260], caminho[260];
	int status, num;
	static FILE *fptr;

	
	DisplayPanel (medicao_painel);
	while (TRUE) {
 		GetUserEvent (1, &objeto, &item);
   		if (objeto == menu_barra) {
			switch (item) {
				case MENU_CONFIGURACAO:
					novo = FALSE;
					break;
				case MENU_MEDICAO:
					novo = FALSE;
					break;
				case MENU_INCERT:
					novo = FALSE;
					break;
				case MENU_GRAFICOS:
					novo = FALSE;
					break;
				case MENU_SAIR:
					novo = TRUE;
					HidePanel (configuracao_painel);
					HidePanel (medicao_painel);
					HidePanel (graficos_painel);
					HidePanel (incert_painel);
					break;
			}
		return;
		}
		if (objeto == medicao_painel) {
			switch (item) {
				case MEDICAO_MEDIR:
					MessagePopup ("NOTE",
						"CHOOSE THE CONVERGENCE PARAMETERS");
					//sampled_data();
					SetCtrlAttribute (medicao_painel, MEDICAO_FREQ_NLS, ATTR_DIMMED, 0);
					break;
				case MEDICAO_FREQ_NLS:
					/* starting guess */
					GetCtrlVal (medicao_painel, MEDICAO_STARTPT_NLS, &Freq_start);	
					startpt_NLS = 4*asin(1)*Freq_start*Tsamp;
					GetCtrlVal (medicao_painel, MEDICAO_IMAX_NLS, &itermax_NLS); /* max # of iterations	     */
					GetCtrlVal (medicao_painel, MEDICAO_RHO_BEGIN_NLS, &rho_NLS); /* stepsize geometric shrink */
					GetCtrlVal (medicao_painel, MEDICAO_EPSMIN_NLS, &epsilon_NLS); /* ending value of stepsize  */
					f1 (startpt_NLS); 
					ampli_est ();
					new_bi[0] = 0;
					for (j=1; j<ifun+1; j++) {
						new_bi[j] = newbi[j-1];
					}	
					iters = 0;
					steplength_NLS = 1;
					omega2 = startpt_NLS;
					while ((iters < itermax_NLS) && (steplength_NLS > epsilon_NLS)) {
						f2 (omega2, new_bi);
						freq_ampli_est ();
						omega2 += rho_NLS*end_bi[0];
						f1 (omega2); 
						ampli_est ();
	   					new_bi[0] = end_bi[0];
	   					for (j=1; j<ifun+1; j++) {
	   						new_bi[j] = newbi[j-1];
	   					}	
	   					iters++;
						SetCtrlVal (medicao_painel, MEDICAO_ITERS, iters);
	   					steplength_NLS = fabs(end_bi[0]);
					}	
   					for (j=1; j<ifun+1; j++) {
  						end_bi[j] = new_bi[j];
   					}	
					Freq_NLS = omega2;
					result_freq_NLS ();
					SetCtrlAttribute (medicao_painel, MEDICAO_AMPL_NLS, ATTR_DIMMED, 0);
					break;
				case MEDICAO_AMPL_NLS:
					ampli_NLS ();
					result_ampl_NLS ();
					SetCtrlAttribute (medicao_painel, MEDICAO_MEDIR_FLICKER, ATTR_DIMMED, 0);
					break;
				case MEDICAO_MEDIR_FLICKER:
					MessagePopup ("NOTE",
						"CHOOSE THE CONVERGENCE PARAMETERS");
					//sampled_data1();
					GenerateData (shape);
					SetCtrlAttribute (medicao_painel, MEDICAO_FREQ, ATTR_DIMMED, 0);
					SetCtrlAttribute (medicao_painel, MEDICAO_FREQ_NLS, ATTR_DIMMED, 1);
					SetCtrlAttribute (medicao_painel, MEDICAO_AMPL_NLS, ATTR_DIMMED, 1);
					break;
				case MEDICAO_FREQ:
					GetCtrlVal (medicao_painel, MEDICAO_NVARS, &nvars);
					 /* starting guess for Student prob function */
					GetCtrlVal (medicao_painel, MEDICAO_STARTPT_0, &Freq_start_carr);	
					GetCtrlVal (medicao_painel, MEDICAO_STARTPT_1, &Fase_start);
					GetCtrlVal (medicao_painel, MEDICAO_STARTPT_2, &DeltaF_start);	
					startpt[0] = 4*asin(1)*Freq_start_carr*Tsamp;
					startpt[1] = Fase_start;
					startpt[2] = 4*asin(1)*DeltaF_start*Tsamp;
					//GetValueFromIndex (medicao_painel, MEDICAO_LIST_FREQ, 0, &omega2);
					GetCtrlVal (medicao_painel, MEDICAO_IMAX, &itermax); /* max # of iterations	     */
					GetCtrlVal (medicao_painel, MEDICAO_RHO_BEGIN, &rho); /* stepsize geometric shrink */
					GetCtrlVal (medicao_painel, MEDICAO_EPSMIN, &epsilon); /* ending value of stepsize  */
	   				jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
	  			    for (k=0; k<nvars; k++) {
	  			    	omega[k] = endpt[k];
	  			    } 	
	  			    result_freq ();
					SetCtrlAttribute (medicao_painel, MEDICAO_AMPL, ATTR_DIMMED, 0);
					break;
				case MEDICAO_AMPL:
					ampli();
					result_ampl ();
					SetCtrlAttribute (medicao_painel, MEDICAO_ARQUIVAR, ATTR_DIMMED, 0);
					SetCtrlAttribute (incert_painel, INCERT_LER, ATTR_DIMMED, 0);
					break;
				case MEDICAO_ARQUIVAR:
		   			GetProjectDir (diretorio);
   					status = FileSelectPopup (diretorio, "datafile.txt", "Datafiles (*.txt)", 
      														"DataFile Storage", VAL_SAVE_BUTTON, 0, 1, 1, 1, caminho);
	           		if (status != VAL_NO_FILE_SELECTED) {
       					fptr = fopen (caminho, "w+");
       					for (k=0; k<(ifun+1)/2; k++)
							num = fprintf (fptr, "%10u\t%10f\t%10f\n", k,Mod_NLS[k]);
							fclose (fptr);  
						}
      		      	break;
			}
		}
	}
	
}


void GenerateData (int shape)
{

	int i, j, k;
	double Gaussrnd[N];
	double Gaussrnd_pha[N];


	switch (shape) {
		case 0 :        /* Sinewave */
		  	for (i=0; i<ino; i++) {
	  	 		data1[i] = 0;
	  	 		Onda_Seno[i] = 0;
	  		}
			GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		GaussNoise (ino, 1, 0, Gaussrnd_pha);
		  	for (i=0; i<ino; i++) {	  
	  	 		Onda_Seno[i] += (Ampl_Gen_Mod/(Freq_Gen_Mod))*sin(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]);
	     		//Onda_Seno[i] += dc_Gen_Mod;
	  		}   
	  		for (i=0; i<ino; i++) {
	     		data1[i] = sqrt(2)*Ampl_Gen*cos((4*asin(1)*Freq_Gen*Tsamp+4*asin(1)*(1/Ampl_Gen_Mod)*Delta_F*Tsamp*dc_Gen_Mod)*i+Fase_Gen+(1/Ampl_Gen_Mod)*Delta_F*(Onda_Seno[i]+Gaussrnd[i]))+Gaussrnd[i];
		 		data1[i] += dc_Gen;
	  		}
			break;
	    case 1 :        /* Square wave */		
	  		for (i=0; i<ino; i++) {
	  	 		data1[i] = 0;
	  	 		Onda_Quad_Gen[i] = 0;
	  		}
	  		GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		GaussNoise (ino, 1, 0, Gaussrnd_pha);
	  		for (i=0; i<ino; i++) {	  
	  	 		for (j=0; j<Nharm_Gen_Mod/2; j++) {	   // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    		k = 2*j+1;
	  	    		Onda_Quad_Gen[i] -= (Ampl_Gen_Mod/(k*Freq_Gen_Mod))*cos(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]))/k;
	  	 		}
	     		Onda_Quad_Gen[i] *= 2/asin(1);
	     		//Onda_Quad_Gen[i] += dc_Gen_Mod;
	  		}   
	  		for (i=0; i<ino; i++) {
	     		data1[i] = sqrt(2)*Ampl_Gen*cos((4*asin(1)*Freq_Gen*Tsamp+4*asin(1)*(asin(1)/(2*Ampl_Gen_Mod))*Delta_F*Tsamp*dc_Gen_Mod)*i+Fase_Gen+(asin(1)/(2*Ampl_Gen_Mod))*Delta_F*(Onda_Quad_Gen[i]+Gaussrnd[i]))+Gaussrnd[i];
		 		data1[i] += dc_Gen;
	  		}
			break;
        case 2 :        /* Triangle wave */
	  		for (i=0; i<ino; i++) {
	  	 		data1[i] = 0;
	  	 		Onda_Triang_Gen[i] = 0;
	  		}
	  		GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  		GaussNoise (ino, 1, 0, Gaussrnd_pha);
	  		for (i=0; i<ino; i++) {	  
	  	 		for (j=0; j<Nharm_Gen_Mod/2; j++) {	   // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    		k = 2*j+1;
	  	    		Onda_Triang_Gen[i] += (Ampl_Gen_Mod/(k*Freq_Gen_Mod))*sin(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]))/pow(k, 2);
	  	    	}
	     		Onda_Triang_Gen[i] *= 2/pow(asin(1), 2);
	     		//Onda_Triang_Gen[i] += dc_Gen_Mod;
	     	}   
	  		for (i=0; i<ino; i++) {
	     		data1[i] = sqrt(2)*Ampl_Gen*cos((4*asin(1)*Freq_Gen*Tsamp+4*asin(1)*(pow(asin(1), 2)/(2*Ampl_Gen_Mod))*Delta_F*Tsamp*dc_Gen_Mod)*i+Fase_Gen+(pow(asin(1), 2)/(2*Ampl_Gen_Mod))*Delta_F*(Onda_Triang_Gen[i]+Gaussrnd[i]))+Gaussrnd[i];
		 		data1[i] += dc_Gen;
	  		}
			break;
	}
	
}
        

//********************MINIMIZAÇÃO DE HOOKE AND JEEVES*******************************************************************************


int hooke(int nvars, double *startpt, double *endpt, double rho, double epsilon, int itermax)

{
	   double	   delta[VARS];
	   double	   newf, fbefore, steplength, tmp;
	   double	   xbefore[VARS], newx[VARS];
	   int		   i, j, keep;
	   int		   iters, iadj;
	   for (i = 0; i < nvars; i++) {
		   newx[i] = xbefore[i] = startpt[i]; 
		   delta[i] = fabs(startpt[i] * rho);
		   if (delta[i] == 0.0)
			   delta[i] = rho;
	   }
	   iadj = 0;
	   steplength = rho;
	   iters = 0;
	   fbefore = f(newx, nvars);
	   newf = fbefore;
	   while ((iters < itermax) && (steplength > epsilon)) {
		   iters++;
		   iadj++;
		   printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);
		   for (j = 0; j < nvars; j++)
			   printf("   x[%2d] = %.4le\n", j, xbefore[j]);
		   /* find best new point, one coord at a time */
		   for (i = 0; i < nvars; i++) {
			   newx[i] = xbefore[i];
		   }
		   newf = best_nearby(delta, newx, fbefore, nvars);
		   /* if we made some improvements, pursue that direction */
		   keep = 1;
		   while ((newf < fbefore) && (keep == 1)) {
			   iadj = 0;
			   for (i = 0; i < nvars; i++) {
				   /* firstly, arrange the sign of delta[] */
				   if (newx[i] <= xbefore[i])
					   delta[i] = 0.0 - fabs(delta[i]);
				   else
					   delta[i] = fabs(delta[i]);
				   /* now, move further in this direction */
				   tmp = xbefore[i];
				   xbefore[i] = newx[i];
				   newx[i] = newx[i] + newx[i] - tmp;
			   }
			   fbefore = newf;
			   newf = best_nearby(delta, newx, fbefore, nvars);
			   /* if the further (optimistic) move was bad.... */
			   if (newf >= fbefore)
				   break;
			   /* make sure that the differences between the new */
			   /* and the old points are due to actual */
			   /* displacements; beware of roundoff errors that */
			   /* might cause newf < fbefore */
			   keep = 0;
			   for (i = 0; i < nvars; i++) {
				   keep = 1;
				   if (fabs(newx[i] - xbefore[i]) >
				       (0.5 * fabs(delta[i])))
					   break;
				   else
					   keep = 0;
			   }
		   }
		   if ((steplength >= epsilon) && (newf >= fbefore)) {
			   steplength = steplength * rho;
			   for (i = 0; i < nvars; i++) {
				   delta[i] *= rho;
			   }
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   endpt[i] = xbefore[i];
	   return (iters);
}


/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double *delta, double *point, double prevbest, int nvars)
	   
{
	   double	   z[VARS];
	   double	   minf, ftmp;
	   int		   i;
	   minf = prevbest;
	   for (i = 0; i < nvars; i++)
		   z[i] = point[i];
	   for (i = 0; i < nvars; i++) {
		   z[i] = point[i] + delta[i];
		   ftmp = f(z, nvars);
		   if (ftmp < minf)
			   minf = ftmp;
		   else {
			   delta[i] = 0.0 - delta[i];
			   z[i] = point[i] + delta[i];
			   ftmp = f(z, nvars);
			   if (ftmp < minf)
				   minf = ftmp;
			   else
				   z[i] = point[i];
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   point[i] = z[i];
	   return (minf);
}


/* function - one frequency */


void f1 (double x_NLS)

{

	   omega2 = x_NLS;

	   if (omega2 == 0) {
	   	   omega2 += 1E-6;
	   }

   	   setGIJ_NLS(omega2);

}


void f2 (double x_NLS, double *bi_NLS)

{
	   int k;
	   
	   omega2 = x_NLS;

	   if (omega2 == 0) {
	   	   omega2 += 1E-6;
	   }

   	   a_NLS[0] = bi_NLS[1];
	   b_NLS[0] = 0;
	   for (k=1; k<ifun; k++) {
	   	  if (k < (ifun+1)/2) {
	  	 	  a_NLS[k] = bi_NLS[k+1];
	  	  }
	  	  else {
	  	 	  b_NLS[k-(ifun-1)/2] = bi_NLS[k+1];
	  	  }
	   }
   	   setGIJC(omega2, a_NLS, b_NLS);

}


/* Student prob function - One frequency */


double f(double *x, int n)

{

   	   int i, j;
   	   
   	   funevals++;

	   nfun = NFUN;

	   for (i=0; i<nvars; i++) {
	   	   omega[i] = x[i];
	   }

	   zloge = 0;
	   setGIJ(omega);
   	   prob();
   	   nstloge = -stloge;
   	   nst = -st;		
	   return nstloge;	   
}



/*void sampled_data(void)

									  
{
	  int i, j, k;
	  double Gaussrnd[N];
	  
	  for (i=0; i<ino; i++) {
	  	 data[i] = 0;
	  	 Onda_Seno[i] = 0;
	  	 //Onda_Triang_Gen[i] = 0;
	  	 //Onda_Quad_Gen[i] = 0;
	  }
	  GetCtrlVal (configuracao_painel, CONFIG_SIGMA, &n_sigma);	  
	  GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  for (i=0; i<ino; i++) {	  
	  	 Onda_Seno[i] = Ampl_Gen_Mod*cos(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod);
	  	 /*for (j=0; j<4*Nharm_Gen_Mod; j++) {	  // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    k = 2*j+1;
	  	    Onda_Triang_Gen[i] += Ampl_Gen_Mod*cos(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod))/pow(k, 2);
	  	    Onda_Quad_Gen[i] += Ampl_Gen_Mod*sin(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp-asin(1)+Fase_Gen_Mod))/k;
	  	 }*/
	     //Onda_Triang_Gen[i] *= 2/pow(asin(1), 2);
	     //Onda_Quad_Gen[i] *= 2/asin(1);
/*	     Onda_Seno[i] += dc_Gen_Mod;
	     //Onda_Triang_Gen[i] += dc_Gen_Mod;
	     //Onda_Quad_Gen[i] += dc_Gen_Mod;
	  }   
	  for (i=0; i<ino; i++) {
		 data[i] = Onda_Seno[i]+Gaussrnd[i];
		 //data[i] = Onda_Triang_Gen[i]+Gaussrnd[i];
		 //data[i] = Onda_Quad_Gen[i]+Gaussrnd[i];
	  }
}*/     



/*void sampled_data1(void)

									  
{
	  int i;
	  double Delta_V, Gaussrnd[N];
	  double Gaussrnd_pha[N];

	  for (i=0; i<ino; i++) {
	  	 data[i] = 0;
	  	 Onda_Seno[i] = 0;
	  	 //Onda_Triang_Gen[i] = 0;
	  	 //Onda_Quad_Gen[i] = 0;
	  }
	  GaussNoise (ino, n_sigma, 0, Gaussrnd);
	  GaussNoise (ino, 1, 0, Gaussrnd_pha);
      GetCtrlVal (configuracao_painel, CONFIG_DELTAF, &Delta_F);		// DESVIO DE FREQUÊNCIA GERADO
	  for (i=0; i<ino; i++) {	  
	  	 Onda_Seno[i] += (Ampl_Gen_Mod/(Freq_Gen_Mod))*sin(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]);
	  	 /*for (j=0; j<4*Nharm_Gen_Mod; j++) {	   // PARA RECONSTRUÇÃO PERFEITA DA ONDA, USAR Nharm_Gen_Mod/2
	  	    k = 2*j+1;
	  	    //Onda_Triang_Gen[i] += (Ampl_Gen_Mod/(k*Freq_Gen_Mod))*sin(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]))/pow(k, 2);
	  	    Onda_Quad_Gen[i] -= (Ampl_Gen_Mod/(k*Freq_Gen_Mod))*cos(k*(4*asin(1)*Freq_Gen_Mod*i*Tsamp+Fase_Gen_Mod+3*Gaussrnd_pha[1]))/k;
	  	 }*/
	     //Onda_Triang_Gen[i] *= 2/pow(asin(1), 2);
	     //Onda_Quad_Gen[i] *= 2/asin(1);
/*	     Onda_Seno[i] += dc_Gen_Mod;
	     //Onda_Triang_Gen[i] += dc_Gen_Mod;
	     //Onda_Quad_Gen[i] += dc_Gen_Mod;
	  }   

	  for (i=0; i<ino; i++) {
	     data[i] = sqrt(2)*Ampl_Gen*cos(4*asin(1)*Freq_Gen*i*Tsamp+Fase_Gen+Delta_F*(Onda_Seno[i]+Gaussrnd[i]))+Gaussrnd[i];
	     //data[i] = sqrt(2)*Ampl_Gen*cos(4*asin(1)*Freq_Gen*i*Tsamp+Fase_Gen+Delta_F*(Onda_Triang_Gen[i]+Gaussrnd[i]))+Gaussrnd[i];
	     //data[i] = sqrt(2)*Ampl_Gen*cos(4*asin(1)*Freq_Gen*i*Tsamp+Fase_Gen+Delta_F*(Onda_Quad_Gen[i]+Gaussrnd[i]))+Gaussrnd[i];
		 data[i] += dc_Gen;
	  }
}*/     



//********************MATRIZ DO PROJETO (GIJ)************************************************************************************


void setGIJ_NLS (double omega_2)

{
	  int i, j;
	  
	  
	  for (i=0; i<ino; i++) {
	  	 GIJ_NLS[i][0] = 1;
	  	 for (j=1; j<=ifun/2; j++) {
	  	 	GIJ_NLS[i][j] = cos(j*omega_2*i);
	  	 	GIJ_NLS[i][ifun/2+j] = sin(j*omega_2*i);
	  	 }
	  }
}


void setGIJC (double omega_2, double *a_, double *b_)

{
	  int i, j;
	  
	  for (i=0; i<ino; i++) {
	     for (j=0; j<=ifun/2+1; j++) {
	     	GIJC[i][j] = 0;
	     }	
	  }
	  for (i=0; i<ino; i++) {
	  	 GIJC[i][1] = 1;
	  	 for (j=1; j<=ifun/2; j++) {
	  	 	GIJC[i][0] += -a_[j]*j*i*sin(j*omega_2*i)+b_[j]*j*i*cos(j*omega_2*i);
	  	 }
	  	 for (j=2; j<=ifun/2+1; j++) {
	  	 	GIJC[i][j] = cos((j-1)*omega_2*i);
	  	 	GIJC[i][ifun/2+j] = sin((j-1)*omega_2*i);
	  	 }
	  }
}



void setGIJ(double *omega)

{
	  int i, k;
	  double Onda[N];
	  
	  for (i=0; i<ino; i++) {
	  	 Onda[i] = 0;
	  }
	  for (i=0; i<ino; i++) {
	  	 for (k=1; k<(ifun+1)/2; k++) {
	     	Onda[i] += (Mod_NLS[k]/(k*Freq_NLS))*sin(k*(Freq_NLS*i+omega[1])+Fase_NLS[k]);
	     }
	     //Onda[i] += Mod_NLS[0];
	  }
	  for (i=0; i<ino; i++) {
	  	 GIJ[i][0] = 1;
 	  	 GIJ[i][1] = cos((omega[0]+(1/Mod_NLS[1])*omega[2]*Mod_NLS[0])*i+(1/Mod_NLS[1])*omega[2]*Onda[i]);
	  	 GIJ[i][2] = sin((omega[0]+(1/Mod_NLS[1])*omega[2]*Mod_NLS[0])*i+(1/Mod_NLS[1])*omega[2]*Onda[i]);
	  }

}


//********************CÁLCULO DA PROBABILIDADE DA FREQUÊNCIA ************************************************************************


void prob(void)

{

      double h1, h2, y2, qq, ahold, y2bar, dif;
      int i, j;
      
	  
	  ortho();
	  
	  h2 = 0;
	  for (j=0; j<nfun; j++) {
	  	 h1 = 0;
	  	 for (i=0; i<ino; i++) {
	  	 	 h1 += data1[i]*HIJ[i][j];
	  	 }
	  	 hi[j] = h1;
	  	 h2 += h1*h1;
	  }					
	  h2bar = h2/nfun;
	  
	  y2 =0;
	  for (i=0;i<ino; i++) {
	  	 y2 += data1[i]*data1[i];
	  }
	  y2bar = y2/ino;
	  
	  qq = 1-h2/y2;
	  if (qq <= 0) {
	  	 qq = 1E-16;
	  }
	  stloge = log(qq)*(nfun-ino)/2;
	  ahold = stloge-zloge;
	  st = 0;
	  if (zloge == 0) {
	  	 st = 0;
	  }									  
	  else {
	  	 st = exp(ahold);
	  }
  	  dif = y2bar-nfun*h2bar/ino;
	  sigma = sqrt(fabs(dif)*ino/(ino-nfun-2));
	  phat = nfun*(sigma*sigma+h2bar)*st;
	  
}


//********************ORTOGONALIZAÇÃO DA MATRIZ DO PROJETO (GIJ)***************************************************************************


void ortho(void)

{
	  double TIJ[NFUN][N], M[NFUN][NFUN], V[NFUN][NFUN], D_norm[NFUN][NFUN], A[N][NFUN]; 
	  double D_vec[NFUN], SqrSumCol[NFUN], norm[NFUN];
	  int i, j;
	  
	  Transpose (GIJ, ino, nfun, TIJ);
	  MatrixMul (TIJ, GIJ, nfun, ino, nfun, M);
	  SymEigenValueVector (M, nfun, 1, D_vector, V);
	  for (j=0; j<nfun; j++) {
	  	 SqrSumCol[j] = 0;
	  }
	  for (j=0; j<nfun; j++) {
	  	 for (i=0; i<nfun; i++) {
	  	 	 SqrSumCol[j] += pow(V[i][j], 2);
	  	 }
	  }
	  for (j=0; j<nfun; j++) {
	  	 norm[j] = 0;
	  }
	  for (j=0; j<nfun; j++) {
	  	 norm[j] = sqrt(SqrSumCol[j]);
	  }
	  for (j=0; j<nfun; j++) {
	  	 for (i=0; i<nfun; i++) {
	  	 	 V_norm[i][j] = V[i][j]/norm[j];
	  	 }														  
	  }
	  for (i=0; i<nfun; i++) {
	  	 D_vec[i] = 0;
	  }
	  for (i=0; i<nfun; i++) {
	     //if (D_vector[i] >= 0) {
	     	 D_vec[i] = sqrt(fabs(D_vector[i]));
	     //}
	  }
	  SpecialMatrix (1, nfun, D_vec, nfun, 0, 0, D_norm);
	  InvMatrix (D_norm, nfun, invD_norm);
	  MatrixMul (GIJ, V_norm, ino, nfun, nfun, A);
	  MatrixMul (A, invD_norm, ino, nfun, nfun, HIJ);
}



void ampli_est (void)

{
      double h1;
  	  int i, j;
  	  
  	  for (i=0; i<ino; i++) {
	  	 for (j=0; j<ifun; j++) {
	  	 	GIJ_array[i*ifun+j] = GIJ_NLS[i][j];
	  	 }
	  }
	  Transpose (GIJ_array, ino, ifun, TIJ_array);
	  MatrixMul (TIJ_array, GIJ_array, ifun, ino, ifun, M_array);
	  InvMatrix (M_array, ifun, invM_array);

	  for (j=0; j<ifun; j++) {
	  	 h1 = 0;
	  	 for (i=0; i<ino; i++) {
	  	 	 h1 += data[i]*GIJ_NLS[i][j];
	  	 }
	  	 hi_NLS[j] = h1;
	  }					
	  MatrixMul (invM_array, hi_NLS, ifun, ifun, 1, newbi);

}


void freq_ampli_est (void)

{	      
      double h1;
  	  int i, j;
		  

      for (i=0; i<ino; i++) {
  	  	 for (j=0; j<ifun+1; j++) {
  	 		GIJC_array[i*(ifun+1)+j] = GIJC[i][j];
  	 	 }
  	  }
  	  Transpose (GIJC_array, ino, ifun+1, TIJC_array);
  	  MatrixMul (TIJC_array, GIJC_array, ifun+1, ino, ifun+1, MC_array);
  	  InvMatrix (MC_array, ifun+1, invMC_array);

      for (j=0; j<ifun+1; j++) {
  	     h1 = 0;
  	  	 for (i=0; i<ino; i++) {
  	 	 	h1 += data[i]*GIJC[i][j];
  	  	 }
  	  	 hiC[j] = h1;
      }					
  	  MatrixMul (invMC_array, hiC, ifun+1, ifun+1, 1, end_bi);
}	   



//********************AMPLITUDES E FASES DOS HARMÔNICOS***********************************************************************************


void ampli_NLS (void)

{

	  int i, j, k;
	  double Fase_ini;
	  double *a_corr_NLS, *b_corr_NLS;	  

	  a_corr_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));	  
	  b_corr_NLS = (double *) malloc (((ifun+1)/2)*sizeof(double));	  

   	  a_NLS[0] = end_bi[1];
	  b_NLS[0] = 0;
	  for (k=1; k<ifun; k++) {
	  	 if (k < (ifun+1)/2) {
	  		a_NLS[k] = end_bi[k+1];
	  	 }
	  	 else {
 	  	    b_NLS[k-(ifun-1)/2] = end_bi[k+1];
	  	 }
	  }
	  Mod_NLS[0] = a_NLS[0];
	  Fase_NLS[0] = atan2 (-b_NLS[0], a_NLS[0]);
	  for (k=1; k<(ifun+1)/2; k++) {
	  	 Mod_NLS[k] = sqrt (pow (a_NLS[k], 2)+pow (b_NLS[k], 2));
	  	 Fase_NLS[k] = atan2 (-b_NLS[k], a_NLS[k]);
	  }
	  for (k=1; k<(ifun+1)/2; k++) {
		 Mod_corr_NLS[k] = 0;
		 Fase_corr_NLS[k] = 0;
	  }
	  Acrms_NLS = 0;
	  for (k=1; k<(ifun+1)/2; k++) {
		 Fase_corr_NLS[k] = 0;
		 a_NLS[k] /= sqrt(2); 
		 b_NLS[k] /= sqrt(2);
		 a_corr_NLS[k] = a_NLS[k]*cos(Fase_corr_NLS[k])+b_NLS[k]*sin(Fase_corr_NLS[k]);   /* CORRIGE A FASE */
		 b_corr_NLS[k] = b_NLS[k]*cos(Fase_corr_NLS[k])-a_NLS[k]*sin(Fase_corr_NLS[k]);
		 Acrms_NLS += pow (a_corr_NLS[k], 2)+pow (b_corr_NLS[k], 2);
	  }
	  Mod_Fundamental_NLS = sqrt (pow (a_corr_NLS[1], 2)+pow (b_corr_NLS[1], 2));
	  Fase_ini = Fase_NLS[1];
	  Fase_ini *= 180/(2*asin(1));		// VALOR DA FASE INICIAL
	  for (k=1; k<(ifun+1)/2; k++) {	 
		 Mod_corr_NLS[k] = 100*sqrt (pow (a_corr_NLS[k], 2)+pow (b_corr_NLS[k], 2))/Mod_Fundamental_NLS;	 // VALOR CORRIGIDO DAS AMPLITUDES 
		 Fase_corr_NLS[k] = Fase_NLS[k];  
		 Fase_corr_NLS[k] *= 180/(2*asin(1));		// VALOR CORRIGIDO DAS FASES
	   	 Fase_corr_NLS[k] += k*(-90-Fase_ini);	/* CORREÇÃO PARA DESLOCAMENTO NO TEMPO (FASE LINEAR), PARA CADA TREM DE IMPULSOS */
		 if (abs(Fase_corr_NLS[k])>360) {
			Fase_corr_NLS[k] = fmod (Fase_corr_NLS[k], 360);  /* REDUÇÃO À PRIMEIRA VOLTA */
		 }
	  }
	  Mod_corr_NLS[0] = 100*Mod_NLS[0]/Mod_Fundamental_NLS;
	  Fase_corr_NLS[0] = 0;
	  Acrms_NLS = sqrt (Acrms_NLS);
	  if (pow(Acrms_NLS, 2) > pow (Mod_Fundamental_NLS, 2)) {
		 THD = 100*(sqrt (pow (Acrms_NLS, 2)-pow (Mod_Fundamental_NLS, 2)))/Mod_Fundamental_NLS;	 
		 															// PROTEÇÃO RAIZ QUADRADA DE NÚMERO NEGATIVO
	  }
	  else {
		 THD = 0;
	  }
	  for (i=0; i<ino; i++) {
	  	 Onda_Result[i] = 0;
	  	 residue[i] = 0;
	  }
	  for (i=0; i<ino; i++) {
	  	 Onda_Result[i] = Mod_NLS[0];
	  	 for (j=1; j<=ifun/2; j++) {
	  	 	Onda_Result[i] += Mod_NLS[j]*cos(j*Freq_NLS*i+Fase_NLS[j]);
	  	 }
	  }
	  for (i=0; i<ino; i++) {
	  	 residue[i] = Onda_Result[i]-data[i];
	  }
  	  for (i=0; i<ino; i++) {
	  	 data_print[i] = data[i];
	  }

}
	  
	  
void ampli (void)

{

	  int i, k;
	  double bi[NFUN];
	  double a[(NFUN+1)/2], b[(NFUN+1)/2];
	  double a_corr[(NFUN+1)/2], b_corr[(NFUN+1)/2];
	  double B[NFUN][NFUN];
	  double  Onda_Result_Temp[N];
	  
	  MatrixMul (V_norm, invD_norm, nfun, nfun, nfun, B);
	  MatrixMul (B, hi, nfun, nfun, 1, bi);
	  a[0] = bi[0];
	  b[0] = 0;
	  for (k=1; k<nfun; k++) {
	  	 if (k < (nfun+1)/2) {
	  	 	 a[k] = bi[k];
	  	 }
	  	 else {
	  	 	 b[k-(nfun-1)/2] = bi[k];
	  	 }
	  }
  	  Modulo[0] = a[0];
  	  Fase[0] = 0;
	  for (k=1; k<(nfun+1)/2; k++) {
		 Modulo[k] = 0;
		 Fase[k] = 0;
	  }
  	  for (k=1; k<(nfun+1)/2; k++) {
  	  	 Modulo[k] = sqrt (pow (a[k], 2)+pow (b[k], 2));	 // VALOR DAS AMPLITUDES
  	  	 Fase[k] = atan2 (-b[k],a[k]);  					 // VALOR DAS FASES
  	  }
	  for (k=1; k<(nfun+1)/2; k++) {
		 a_corr[k] = 0;
		 b_corr[k] = 0;
		 Modulo_corr[k] = 0;
		 Fase_corr[k] = 0;
	  }
  	  Acrms = 0;
  	  for (k=1; k<(nfun+1)/2; k++) {
		 a[k] /= sqrt (2);
		 b[k] /= sqrt (2);
		 Fase_corr[k] = 0;
		 a_corr[k] = a[k]*cos(Fase_corr[k])+b[k]*sin(Fase_corr[k]);   /* CORRIGE A FASE */	 
		 b_corr[k] = b[k]*cos(Fase_corr[k])-a[k]*sin(Fase_corr[k]);
		 Acrms += pow (a_corr[k], 2)+pow (b_corr[k], 2);
	  }
	  Mod_Fundamental = sqrt (pow (a_corr[1], 2)+pow (b_corr[1], 2));
	  Acrms = sqrt (Acrms);
	  for (k=1; k<(nfun+1)/2; k++) {	 
		 Modulo_corr[k] = 100*sqrt (pow (a_corr[k], 2)+pow (b_corr[k], 2))/Mod_Fundamental;	 // VALOR CORRIGIDO DAS AMPLITUDES
		 Fase_corr[k] = Fase[k];  
		 Fase_corr[k] *= 180/(2*asin(1));		// VALOR DAS FASES
	  }
	  Modulo_corr[0] = 100*Modulo[0]/Mod_Fundamental;
	  Fase_corr[0] = 0;
  	  for (i=0; i<ino; i++) {
	  	 Onda_Result[i] = 0;
	  	 Onda_Result_Temp[i] = 0;
	  }
	  for (i=0; i<ino; i++) {
	  	 for (k=1; k<(ifun+1)/2; k++) {
	     	Onda_Result[i] += (Mod_NLS[k]/(k*Freq_NLS))*sin(k*(Freq_NLS*i+omega[1])+Fase_NLS[k]);
	     }
	     //Onda_Result[i] += Mod_NLS[0];
	  }
	  for (i=0; i<ino; i++) {
		 residue[i] = Modulo[0]+Modulo[1]*cos((omega[0]+(1/Mod_NLS[1])*omega[2]*Mod_NLS[0])*i+(1/Mod_NLS[1])*omega[2]*Onda_Result[i]+Fase[1])-data1[i];
	  }
 	  
	  for (i=0; i<ino; i++) {
		 Onda_Result[i] *= (1/Mod_NLS[1])*omega[2];
		 Onda_Result[i] += (omega[0]+(1/Mod_NLS[1])*omega[2]*Mod_NLS[0])*i+Fase[1];
		 Onda_Result_Temp[i] = cos(Onda_Result[i]);
		 Onda_Result[i] = Onda_Result_Temp[i];
		 Onda_Result[i] *= Modulo[1];
		 Onda_Result[i] += Modulo[0];
	  }
	  for (i=0; i<ino; i++) {
	  	 data_print[i] = data1[i];
	  }
	  
}

//********************MOSTRA RESULTADOS DE FREQUÊNCIA NLS ****************************************************************************


void result_freq_NLS (void)

{
	  int k;
	  char harm_freq[50];

	  SetCtrlVal (medicao_painel, MEDICAO_ENDPT_NLS, Freq_NLS/(4*asin(1)*Tsamp));
	  ClearListCtrl (medicao_painel, MEDICAO_LIST_FREQ);
	  for (k=0; k<=ifun/2; k++) {
	  	 Fmt (harm_freq, "%s<%d \033p30 %f[p13]", k, k*Freq_NLS);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_FREQ, -1, harm_freq, k*Freq_NLS);
	  }

}


//********************MOSTRA RESULTADOS DE FREQUÊNCIA ****************************************************************************


void result_freq (void)

{
	  int k;
	  char harm_freq[50];

	  SetCtrlVal (medicao_painel, MEDICAO_ITERS, jj);
	  SetCtrlVal (medicao_painel, MEDICAO_NSTLOGE, nstloge);
	  SetCtrlVal (medicao_painel, MEDICAO_ENDPT_0, omega[0]/(4*asin(1)*Tsamp));
	  SetCtrlVal (medicao_painel, MEDICAO_ENDPT_1, omega[1]);  
	  SetCtrlVal (medicao_painel, MEDICAO_ENDPT_2, omega[2]/(4*asin(1)*Tsamp));
	  ClearListCtrl (medicao_painel, MEDICAO_LIST_FREQ);
	  for (k=0; k<nvars; k++) {
	  	 Fmt (harm_freq, "%s<%d \033p30 %f[p11]", k, omega[k]);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_FREQ, -1, harm_freq, omega[k]);
	  }

}


//********************MOSTRA RESULTADOS DE AMPLITUDES E FASES DE HARMÔNICOS NLS ********************************************************


void result_ampl_NLS (void)

{
	  int k;
	  char harm_amp[50], harm_fase[50];

	  ClearListCtrl (medicao_painel, MEDICAO_LIST_AMP);
	  ClearListCtrl (medicao_painel, MEDICAO_LIST_FASE); 
	  for (k=0; k<(ifun+1)/2; k++) {
	  	 Fmt (harm_amp, "%s<%d \033p30 %f[p5]", k, Mod_corr_NLS[k]);
	  	 Fmt (harm_fase, "%s<%d \033p30 %f[p4]", k, Fase_corr_NLS[k]);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_AMP, -1, harm_amp, Mod_corr_NLS[k]);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_FASE, -1, harm_fase, Fase_corr_NLS[k]);
	  }
	  SetCtrlVal (medicao_painel, MEDICAO_FUNDAMENTAL_NLS, Mod_Fundamental_NLS);
	  SetCtrlVal (medicao_painel, MEDICAO_ACRMS_NLS, Acrms_NLS);
	  SetCtrlVal (medicao_painel, MEDICAO_THD, THD);

}


//********************MOSTRA RESULTADOS DE AMPLITUDES E FASES ********************************************************


void result_ampl (void)

{
	  int k;
	  char harm_amp[50], harm_fase[50];

	  ClearListCtrl (medicao_painel, MEDICAO_LIST_AMP);
	  ClearListCtrl (medicao_painel, MEDICAO_LIST_FASE); 
	  for (k=0; k<(nfun+1)/2; k++) {
	  	 Fmt (harm_amp, "%s<%d \033p30 %f[p5]", k, Modulo_corr[k]);
	  	 Fmt (harm_fase, "%s<%d \033p30 %f[p4]", k, Fase_corr[k]);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_AMP, -1, harm_amp, Modulo_corr[k]);
	  	 InsertListItem (medicao_painel, MEDICAO_LIST_FASE, -1, harm_fase, Fase_corr[k]);
	  }
	  SetCtrlVal (medicao_painel, MEDICAO_FUNDAMENTAL, Mod_Fundamental);
}


//********************PLOTA OS GRÁFICOS ***************************************************************************************


void Plotar (void)

{

	int xmin, xmax, xmin_res, xmax_res;
	int i, j, id, imax, imin, jmax, jmin;
	double ymin, ymax, ymin_res, ymax_res;

	DisplayPanel (graficos_painel);
	while (TRUE) {
 		GetUserEvent (1, &objeto, &item);
   		if (objeto == menu_barra) {
			switch (item) {
				case MENU_CONFIGURACAO:
					novo = FALSE;
					break;
				case MENU_MEDICAO:
					novo = FALSE;
					break;
				case MENU_INCERT:
					novo = FALSE;
					break;
				case MENU_GRAFICOS:
					novo = FALSE;
					break;
				case MENU_SAIR:
					novo = TRUE;
					HidePanel (configuracao_painel);
					HidePanel (medicao_painel);
					HidePanel (graficos_painel);
					HidePanel (incert_painel);
					break;
			}
		return;
		}
		if (objeto == graficos_painel) {
			switch (item) {
				case GRAFICOS_DADOS:
					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_DATA,-1, VAL_IMMEDIATE_DRAW);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_DATA, VAL_XAXIS, VAL_MANUAL, 0, ino);
					PlotY (graficos_painel, GRAFICOS_GRAPH_DATA, data_print, ino, VAL_DOUBLE,
						   VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1, VAL_BLACK);

					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_SIGNAL,-1, VAL_IMMEDIATE_DRAW);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_SIGNAL, VAL_XAXIS, VAL_MANUAL, 0, ino);
					PlotY (graficos_painel, GRAFICOS_GRAPH_SIGNAL, Onda_Result,
						   ino, VAL_DOUBLE, VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1,
						   VAL_BLACK);

					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_RESIDUE,-1, VAL_IMMEDIATE_DRAW);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_RESIDUE, VAL_XAXIS, VAL_MANUAL, 0, ino);
					PlotY (graficos_painel, GRAFICOS_GRAPH_RESIDUE, residue, ino,
						   VAL_DOUBLE, VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1,
						   VAL_BLACK);
					break;
				case GRAFICOS_ESCALA_DADOS:
					GetCtrlVal (graficos_painel, GRAFICOS_XMIN, &xmin);
					GetCtrlVal (graficos_painel, GRAFICOS_XMAX, &xmax);
					GetCtrlVal (graficos_painel, GRAFICOS_YMIN, &ymin);
					GetCtrlVal (graficos_painel, GRAFICOS_YMAX, &ymax);

					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_DATA, VAL_XAXIS, VAL_MANUAL, xmin, xmax);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_DATA, VAL_LEFT_YAXIS, VAL_MANUAL, ymin, ymax);
					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_DATA,-1, VAL_IMMEDIATE_DRAW);					
					PlotY (graficos_painel, GRAFICOS_GRAPH_DATA, data_print, ino, VAL_DOUBLE,
						   VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1, VAL_BLACK);

					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_SIGNAL, VAL_XAXIS, VAL_MANUAL, xmin, xmax);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_SIGNAL, VAL_LEFT_YAXIS, VAL_MANUAL, ymin, ymax);
					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_SIGNAL,-1, VAL_IMMEDIATE_DRAW);					
					PlotY (graficos_painel, GRAFICOS_GRAPH_SIGNAL, Onda_Result,
						   ino, VAL_DOUBLE, VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1,
						   VAL_BLACK);
					
					GetCtrlVal (graficos_painel, GRAFICOS_XMIN_RES, &xmin_res);
					GetCtrlVal (graficos_painel, GRAFICOS_XMAX_RES, &xmax_res);
					GetCtrlVal (graficos_painel, GRAFICOS_YMIN_RES, &ymin_res);
					GetCtrlVal (graficos_painel, GRAFICOS_YMAX_RES, &ymax_res);

					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_RESIDUE, VAL_XAXIS, VAL_MANUAL, xmin_res, xmax_res);
					SetAxisScalingMode (graficos_painel, GRAFICOS_GRAPH_RESIDUE, VAL_LEFT_YAXIS, VAL_MANUAL, ymin_res, ymax_res);
					DeleteGraphPlot (graficos_painel, GRAFICOS_GRAPH_RESIDUE,-1, VAL_IMMEDIATE_DRAW);					
					PlotY (graficos_painel, GRAFICOS_GRAPH_RESIDUE, residue, ino,
						   VAL_DOUBLE, VAL_SCATTER, VAL_SIMPLE_DOT, VAL_SOLID, 1,
						   VAL_BLACK);
					break;
			}
		}

	}
	
}


//********************AVALIA A INCERTEZA DA MEDIÇÃO*******************************************************************************


void Incerteza (void)

{
	int i,j, k, nrows;
	char b;
	double sigma_err, delta;
	double omega_jpkp[3], omega_jmkp[3], omega_jpkm[3], omega_jmkm[3], omega_jp[3], omega_jm[3];
	double H2bar_jpkp[3][3], H2bar_jmkp[3][3], H2bar_jpkm[3][3], H2bar_jmkm[3][3];
	double h2bar_jp[3], h2bar_jm[3];
	double D[NFUN][NFUN], invD[NFUN][NFUN], T_norm[NFUN][NFUN], C[NFUN][NFUN], COV_ampl[NFUN][NFUN]; 
	double Diag[3][3], invDiag[3][3], T_eignorm[3][3], E[3][3], COV_freq[3][3];

	char diretorio[260], caminho[260];
	int status, num;
	static FILE *fptr;
	double Err_Mod_Abs;

	DisplayPanel (incert_painel);
	while (TRUE) {
 		GetUserEvent (1, &objeto, &item);
   		if (objeto == menu_barra) {
			switch (item) {
				case MENU_CONFIGURACAO:
					novo = FALSE;
					break;
				case MENU_MEDICAO:
					novo = FALSE;
					break;
				case MENU_INCERT:
					novo = FALSE;
					break;
				case MENU_GRAFICOS:
					novo = FALSE;
					break;
				case MENU_SAIR:
					novo = TRUE;
					HidePanel (configuracao_painel);
					HidePanel (medicao_painel);
					HidePanel (graficos_painel);
					HidePanel (incert_painel);
					break;
			}
		return;
		}
		if (objeto == incert_painel) {
			switch (item) {
				case INCERT_LER:
					SetCtrlVal (incert_painel, INCERT_LED, TRUE);
					sigma_err = sqrt(pow(sigma, 2)*2/(ino-nfun-4));
				    /*for (i=0; i<nfun; i++) {
	     			   D_vector[i] = fabs(D_vector[i]);
	     			}*/
					SpecialMatrix (1, nfun, D_vector, nfun, 0, 0, D);
	  				InvMatrix (D, nfun, invD);
	  				Transpose (V_norm, nfun, nfun, T_norm);
					MatrixMul (V_norm, invD, nfun, nfun, nfun, C);	  				
					MatrixMul (C, T_norm, nfun, nfun, nfun, COV_ampl);
				  	for (j=0; j<nfun; j++) {
	  	 				for (i=0; i<nfun; i++) {
	  	 	 				COV_ampl[i][j] *= pow(sigma, 2);
	  	 				}														  
	  				}

					delta = 1E-6;
					for (i=0; i<3; i++) {
						omega_jpkp[i] = omega[i];
						omega_jmkp[i] = omega[i];
						omega_jpkm[i] = omega[i];
						omega_jmkm[i] = omega[i];
					}
					for (j=0; j<3; j++) {
						for (k=0; k<3; k++) {
							H2bar_jpkp[j][k] = 0;
							H2bar_jmkp[j][k] = 0;
							H2bar_jpkm[j][k] = 0;
							H2bar_jmkm[j][k] = 0;
							Dh2bar_djdk[j][k] = 0;
						}
					}
					for (j=0; j<3; j++) {
						for (k=0; k<3; k++) {
							if (k>j) {
								omega_jpkp[j] = omega[j]+delta*omega[j];
								omega_jpkp[k] = omega[k]+delta*omega[k];
								setGIJ(omega_jpkp);
								prob ();
								H2bar_jpkp[j][k] = h2bar;
								omega_jmkp[j] = omega[j]-delta*omega[j];
								omega_jmkp[k] = omega[k]+delta*omega[k];
								setGIJ(omega_jmkp);
								prob ();
								H2bar_jmkp[j][k] = h2bar;
								omega_jpkm[j] = omega[j]+delta*omega[j];
								omega_jpkm[k] = omega[k]-delta*omega[k];
								setGIJ(omega_jpkm);
								prob ();
								H2bar_jpkm[j][k] = h2bar;
								omega_jmkm[j] = omega[j]-delta*omega[j];
								omega_jmkm[k] = omega[k]-delta*omega[k];
								setGIJ(omega_jmkm);
								prob ();
								H2bar_jmkm[j][k] = h2bar;
								Dh2bar_djdk[j][k] = (H2bar_jpkp[j][k]-H2bar_jmkp[j][k]-H2bar_jpkm[j][k]+H2bar_jmkm[j][k])
																/(4*pow(delta, 2)*omega[j]*omega[k]);			
								Dh2bar_djdk[k][j] = Dh2bar_djdk[j][k]; 
								for (i=0; i<3; i++) {
									omega_jpkp[i] = omega[i];
									omega_jmkp[i] = omega[i];
									omega_jpkm[i] = omega[i];
									omega_jmkm[i] = omega[i];
								}
							}
						}	
					}		
					for (i=0; i<3; i++) {
						omega_jp[i] = omega[i];
						omega_jm[i] = omega[i];
					}
					for (j=0; j<3; j++) {
						h2bar_jp[j] = 0;
						h2bar_jm[j] = 0;
					}
					for (j=0; j<3; j++) {
						omega_jp[j] = omega[j]+delta*omega[j];
						setGIJ(omega_jp);
						prob ();
						h2bar_jp[j] = h2bar;
						omega_jm[j] = omega[j]-delta*omega[j];
						setGIJ(omega_jm);
						prob ();
						h2bar_jm[j] = h2bar;
						setGIJ(omega);
						prob ();
						Dh2bar_djdk[j][j] = (h2bar_jp[j]-2*h2bar+h2bar_jm[j])/pow(delta*omega[j], 2);
						for (i=0; i<3; i++) {
							omega_jp[i] = omega[i];
							omega_jm[i] = omega[i];
						}
					}		
					for (j=0; j<3; j++) {
						for (k=0; k<3; k++) {
							Dh2bar_djdk[j][k] *= -nfun;
							Dh2bar_djdk[j][k] /= 2;
						}
					}
					/*ortho_err ();
					SpecialMatrix (1, 3, D_eig, 3, 0, 0, Diag);
					InvMatrix (Diag, 3, invDiag);
					Transpose (V_eignorm, 3, 3, T_eignorm);
					MatrixMul (V_eignorm, invDiag, 3, 3, 3, E);
					MatrixMul (E, T_eignorm, 3, 3, 3, COV_freq);*/
				  	
				  	InvMatrix (Dh2bar_djdk, 3, COV_freq);
				  	
				  	for (j=0; j<3; j++) {
	  	 				for (i=0; i<3; i++) {
	  	 	 				COV_freq[i][j] *= pow(sigma, 2);
	  	 				}														  
	  				}
					COV_freq[0][0] /= pow(4*asin(1)*Tsamp, 2);
					COV_freq[2][2] /= pow(4*asin(1)*Tsamp, 2);
					SetCtrlVal (incert_painel, INCERT_SIGMA, sigma);
					SetCtrlVal (incert_painel, INCERT_SIGMA_ERR, sigma_err);
					GetNumTableRows (incert_painel, INCERT_COVAR_AMP, &nrows);
					if (nrows > 0) {
						DeleteTableRows (incert_painel, INCERT_COVAR_AMP, 1, nfun);
						DeleteTableColumns (incert_painel, INCERT_COVAR_AMP, 1, nfun);						
						DeleteTableRows (incert_painel, INCERT_COVAR_FREQ, 1, 3);
						DeleteTableColumns (incert_painel, INCERT_COVAR_FREQ, 1, 3);
					}
					InsertTableRows (incert_painel, INCERT_COVAR_AMP, 1, nfun, VAL_CELL_NUMERIC);
					InsertTableColumns (incert_painel, INCERT_COVAR_AMP, 1, nfun,
									VAL_CELL_NUMERIC);
					InsertTableRows (incert_painel, INCERT_COVAR_FREQ, 1, 3, VAL_CELL_NUMERIC);
					InsertTableColumns (incert_painel, INCERT_COVAR_FREQ, 1, 3,
									VAL_CELL_NUMERIC);
					SetTableRowAttribute (incert_painel, INCERT_COVAR_AMP, -1,
										  ATTR_USE_LABEL_TEXT, 0);
					SetTableColumnAttribute (incert_painel, INCERT_COVAR_AMP, -1,
											 ATTR_USE_LABEL_TEXT, 0);
					SetTableRowAttribute (incert_painel, INCERT_COVAR_FREQ, -1,
										  ATTR_USE_LABEL_TEXT, 0);
					SetTableColumnAttribute (incert_painel, INCERT_COVAR_FREQ, -1,
											 ATTR_USE_LABEL_TEXT, 0);
				  	for (j=1; j<=nfun; j++) {
	  	 				for (i=1; i<=nfun; i++) {
							SetTableCellVal (incert_painel, INCERT_COVAR_AMP, MakePoint (i,j),
											 COV_ampl[i-1][j-1]);
	  	 				}														  
	  				}
				  	for (j=1; j<=3; j++) {
	  	 				for (i=1; i<=3; i++) {
							SetTableCellVal (incert_painel, INCERT_COVAR_FREQ, MakePoint (i,j),
											 COV_freq[i-1][j-1]);
	  	 				}														  
	  				}
		 			SetCtrlVal (incert_painel, INCERT_LED, FALSE);
		 			SetCtrlAttribute (incert_painel, INCERT_ARQUIVAR, ATTR_DIMMED, 0);
					break;
				case INCERT_ARQUIVAR:
	  				GetProjectDir (diretorio);
					status = FileSelectPopup (diretorio, "datafile.txt", "Datafiles (*.txt)", 
														"DataFile Storage", VAL_SAVE_BUTTON, 0, 1, 1, 1, caminho);
					if (status != VAL_NO_FILE_SELECTED) {
   						fptr = fopen (caminho, "w+");
   						for (k=1; k<(nfun+1)/2; k++)
							num = fprintf (fptr, "%10u\t%10f\n", k,Err_Mod_Abs);
							fclose (fptr);  
						}
					break;
			}
		}

	}


}


/*void ortho_err(void)

{
	  double V_eig[3][3]; 
	  double SqrSumCol[NFUN], norm[NFUN];
	  int i, j;
	  
	  SymEigenValueVector (Dh2bar_djdk, 3, 1, D_eig, V_eig);
	  for (j=0; j<3; j++) {
	  	 SqrSumCol[j] = 0;
	  }
	  for (j=0; j<3; j++) {
	  	 for (i=0; i<3; i++) {
	  	 	 SqrSumCol[j] += pow(V_eig[i][j], 2);
	  	 }
	  }
	  for (j=0; j<3; j++) {
	  	 norm[j] = 0;
	  }
	  for (j=0; j<3; j++) {
	  	 norm[j] = sqrt(SqrSumCol[j]);
	  }
	  for (j=0; j<3; j++) {
	  	 for (i=0; i<3; i++) {
	  	 	 V_eignorm[i][j] = V_eig[i][j]/norm[j];
	  	 }														  
	  }
}*/


