

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "gmp.h"



void pgcd_it(mpz_t gcd, mpz_t a, mpz_t b)//version itérative
{
	mpz_t t,copie_a,copie_b;
	mpz_inits(t,copie_a,copie_b,NULL);
	mpz_set(copie_a,a);
	mpz_set(copie_b,b);
	while(mpz_cmp_ui(copie_b,0)!=0){
		mpz_set(t,copie_b);
		mpz_mod(copie_b,copie_a,copie_b);
		mpz_set(copie_a,t);
	}
	mpz_set(gcd,copie_a);
	mpz_clears(t,copie_a,copie_b,NULL);
}

int invert(mpz_t u, mpz_t a, mpz_t n){
	mpz_t gcd;
	mpz_init(gcd);
	pgcd_it(gcd,a,n);
	if(mpz_cmp_ui(gcd,1) != 0)
	{
		printf("Elément non inversible !\n\n");
		mpz_clear(gcd);
		return 0;
	}

	mpz_t u1,v,v1,r,r1,q,rs,us,vs,tmp1,tmp2,tmp3;
	mpz_inits(u1,v,v1,r,r1,q,rs,us,vs,tmp1,tmp2,tmp3,NULL);
	//init
	mpz_set(r,a); // r = a
	mpz_set(r1,n); //r' = n
	mpz_set_ui(u,1); //u = 1
	mpz_set_ui(v,0); //v = 0
	mpz_set_ui(u1,0); //u' = 0
	mpz_set_ui(v1,1); //v' = 1
	//exec
	while(mpz_cmp_ui(r1,0)!=0){
		mpz_fdiv_q(q,r,r1);
		mpz_set(rs,r); // rs = r
		mpz_set(us,u); // us = u
		mpz_set(vs,v); //vs = v
		mpz_set(r,r1); // r = r1
		mpz_set(u,u1); //u = u1
		mpz_set(v,v1); 
		mpz_mul(tmp1,q,r1);
		mpz_mul(tmp2,q,u1);
		mpz_mul(tmp3,q,v1);
		mpz_sub(r1,rs,tmp1);
		mpz_sub(u1,us,tmp2);
		mpz_sub(v1,vs,tmp3);
	}
	mpz_mod(u,u,n);
	mpz_clears(u1,v,v1,r,r1,q,rs,us,vs,tmp1,tmp2,tmp3,NULL);
	return 1;

}

void joye_ladder(mpz_t result, mpz_t m, mpz_t d, mpz_t n){
	int i,k,bit_d;
	mpz_t tab[2];
	mpz_inits(tab[0],tab[1],NULL);
	mpz_set_ui(tab[0],1);
	mpz_set(tab[1],m);
	k = mpz_sizeinbase(d,2);
	for(i = 0; i < k; i++){
		bit_d = mpz_tstbit(d,i);
		mpz_mul(tab[1-bit_d],tab[1-bit_d],tab[1-bit_d]);
		mpz_mul(tab[1-bit_d],tab[1-bit_d],tab[bit_d]);
		mpz_mod(tab[1-bit_d],tab[1-bit_d],n);
	}
	mpz_set(result,tab[0]);
	mpz_clears(tab[0],tab[1],NULL);
}

int test_miller_rabin(mpz_t n, int t){
	if (mpz_cmp_ui(n,0) == 0 || mpz_cmp_ui(n,1) == 0 || mpz_cmp_ui(n,4) == 0){
		return 0;
	}
	mpz_t a,y,n1,n_1,n_2,r,pui,square; //pui = 2^s
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	long int s = 0;
	int j;
	mpz_inits(a,y,n_1,n1,n_2,r,pui,square,NULL);
	mpz_set_ui(pui,1);
	mpz_set_ui(square,2);
	mpz_sub_ui(n_1,n,1);
	mpz_set(n1,n_1);
	mpz_sub_ui(n_2,n_1,1);
	while(!mpz_tstbit(n1,0)){ //tant que c'est pair
		s++;
		mpz_mul_ui(pui,pui,2); //pui = pui*2
		mpz_fdiv_q_2exp(n1,n1,1); //div par 2
	}
 	mpz_divexact(r,n_1,pui);
 	for (int i = 0; i < t; i++){
		while (mpz_cmp_ui(a,0) == 0 || mpz_cmp_ui(a,1) == 0){
			mpz_urandomm(a,generateur,n_2);
		}
		joye_ladder(y,a,r,n);
 		if (mpz_cmp_ui(y,1) != 0 && mpz_cmp(y,n_1) != 0){
 			j = 1;
 			while ((j <= s-1) && (mpz_cmp(y,n_1))){
 				joye_ladder(y,y,square,n);
 				if(mpz_cmp_ui(y,1) == 0){
 					return 0; //composé
 				}
 				j++;
 			}
 			if (mpz_cmp(y,n_1)!= 0){
 				return 0;
 			}
		}
	}
	mpz_clears(a,y,n1,n_1,n_2,r,pui,NULL);
	return (1);
}

void est_egale(mpz_t m,mpz_t m_obtenu){
	if (mpz_cmp (m_obtenu,m) == 0)
	{
		gmp_printf("égalité vérifiée %Zx = %Zx\n\n",m,m_obtenu);
	}
	else
	{
		printf("Erreur\n\n");
	}
}

unsigned long int* k_plus_petit_premier(int k){
	mpz_t z_p,z_next;
	unsigned long int* tab = NULL;
	mpz_inits(z_p,z_next,NULL);
	tab = malloc(sizeof(unsigned long int)*k);
	tab[0] = 2;
	mpz_set_ui(z_p,2);
	for(int i = 1; i < k; i++){
		mpz_nextprime(z_next,z_p);
		tab[i] = mpz_get_ui(z_next);
		mpz_set(z_p,z_next);
	}
	return tab;
	mpz_clears(z_p,z_next,NULL);
	free(tab);
}

int residus_nuls(mpz_t* tab_residus,int k){
	for (int i = 0; i < k-1; i++){
		if (tab_residus[i] == 0){
			return 1;
		}
	}
	return 0;	
}

void crible_opti(mpz_t q, unsigned int b, int k, unsigned int* tab){
	//mpz_set_ui(q,0);
	mpz_t z_residus;
	mpz_t* tab_residus;
	int trouve = 0;
	int div = 0; // = à 1 si au moins un résidu du tableau est nul (donc q divisible et non premier), 0 sinon
	int j;

	mpz_init(z_residus);
	tab_residus = (mpz_t*) malloc((k-1)*sizeof(mpz_t));

	for (j = 0; j< k-1; j++){
		mpz_init(tab_residus[j]);
		mpz_mod_ui(tab_residus[j],q,tab[j+1]);
		if (tab_residus[j] == 0){
			div = 1;
		}
	}

	while(!trouve){
		while(div){
			for(j = 1; j < k ; j++){
				mpz_add_ui(tab_residus[j-1],tab_residus[j-1],2);
				mpz_mod_ui(tab_residus[j-1],tab_residus[j-1],tab[j]);
			}
			mpz_add_ui(q,q,2);
			div = residus_nuls(tab_residus,k);
		}
		if(test_miller_rabin(q,3) == 0){ // si q est quand même composé
			for (j = 1; j < k ; j++){
				mpz_add_ui(tab_residus[j-1],tab_residus[j-1],2);
				mpz_mod_ui(tab_residus[j-1],tab_residus[j-1],tab[j]);
			}
			mpz_add_ui(q,q,2);
		}
		else{
			trouve = 1;
		}
	}
	mpz_clears(z_residus,NULL);
	for(j = 0; j < k-1; j++){
    	mpz_clear(tab_residus[j]);
    }
    free(tab_residus);
}

void generation_premier_valable(mpz_t p,mpz_t phip, gmp_randstate_t generateur, mpz_t e, unsigned int long taille_bit, int k, unsigned int* tab){
	mpz_t gcd;
	mpz_init(gcd);

	do{
		mpz_urandomb(p,generateur,taille_bit);

		mpz_setbit(p,taille_bit - 1); //le forcer a être de taille b
		mpz_setbit(p,0);//le forcer à être impair

    	crible_opti(p,taille_bit,k,tab);

		mpz_sub_ui(phip,p,1);

		pgcd_it(gcd,e,phip);
	}while(mpz_cmp_ui(gcd,1) != 0);
	mpz_clear(gcd);
}

void generation_RSA(mpz_t n, mpz_t e, mpz_t d, int taille_bit, unsigned int* tab, int k){
	mpz_t p,q,phi,phip,phiq;
	mpz_inits(p,q,phi,phip,phiq,NULL);
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	do{
		if(taille_bit%2 == 1){
			generation_premier_valable(p,phip,generateur,e,taille_bit/2,k,tab);
			generation_premier_valable(q,phiq,generateur,e,taille_bit/2 + 1,k,tab);
		}else{
			generation_premier_valable(p,phip,generateur,e,taille_bit/2,k,tab);
			generation_premier_valable(q,phiq,generateur,e,taille_bit/2,k,tab);
		}
		mpz_mul(n,p,q);
	}while(mpz_sizeinbase(n,2)!=taille_bit);
	mpz_mul(phi,phip,phiq);
	mpz_invert(d,e,phi);
	mpz_inits(p,q,phi,phip,phiq,NULL);
	gmp_randclear(generateur);
}

void generation_RSA_CRT(mpz_t n, mpz_t e, mpz_t d_p, mpz_t d_q, mpz_t I_p, mpz_t p, mpz_t q, unsigned long int taille_bit, unsigned int* tab, int k){
	mpz_t phi,phip,phiq;
	mpz_inits(phi,phip,phiq,NULL);
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	do{
		if(taille_bit%2 == 1){
			generation_premier_valable(p,phip,generateur,e,taille_bit/2,k,tab);
			generation_premier_valable(q,phiq,generateur,e,taille_bit/2 + 1,k,tab);
		}else{
			generation_premier_valable(p,phip,generateur,e,taille_bit/2,k,tab);
			generation_premier_valable(q,phiq,generateur,e,taille_bit/2,k,tab);
		}
		mpz_mul(n,p,q);
	}while(mpz_sizeinbase(n,2)!=taille_bit);

	invert(d_p,e,phip);
	invert(d_q,e,phiq);
	invert(I_p,p,q);
	gmp_printf("p = %Zu\nq = %Zu\n\n",p,q);
	gmp_printf("n = 0x%Zx\n\n",n);
	gmp_printf("d_p = 0x%Zx\nd_q = %Zu\n\n",d_p,d_q);
	gmp_printf("I_p = 0x%Zx\n\n",I_p);
	mpz_clears(phip,phiq,phi,NULL);
}


void decrypt_rsa_CRT(mpz_t m, mpz_t c, mpz_t d_p, mpz_t d_q, mpz_t I_p, mpz_t p, mpz_t q){
	mpz_t m_p,m_q;
	mpz_inits(m_p,m_q,NULL);
	joye_ladder(m_p,c,d_p,p);
	joye_ladder(m_q,c,d_q,q);
	mpz_sub(m,m_q,m_p);
	mpz_mul(m,m,I_p);
	mpz_mod(m,m,q);
	mpz_mul(m,p,m);
	mpz_add(m,m_p,m);
	gmp_printf("m_obtenu = 0x%Zx\n\n",m);
	mpz_clears(m_p,m_q,NULL);
}

void sign_CRT(mpz_t sig_CRT, mpz_t m, mpz_t p, mpz_t q, mpz_t d_p, mpz_t d_q,mpz_t I_p,mpz_t n)
{
	mpz_t m_p,m_q,s_p,s_q;
	mpz_inits(m_p,m_q,s_p,s_q,NULL);
	mpz_mod(m_p,m,p);
	mpz_mod(m_q,m,q);
	joye_ladder(s_p,m_p,d_p,p);
	joye_ladder(s_q,m_q,d_q,q);
	mpz_sub(sig_CRT,s_q,s_p);
	mpz_mul(sig_CRT,sig_CRT,I_p);
	mpz_mod(sig_CRT,sig_CRT,q);
	mpz_mul(sig_CRT,sig_CRT,p);
	mpz_add(sig_CRT,sig_CRT,s_p);
	mpz_mod(sig_CRT,sig_CRT,n);
	mpz_clears(m_p,m_q,s_p,s_q,NULL);
}

void messageATraiter(mpz_t m)
{
	int message;
	printf("Entrez le message à chiffrer ou signer\n\n");
	scanf("%d",&message);
	mpz_set_ui(m,message);
	printf("\n");
}



int taille_N()
{
		int Tn_int;
		printf("Donnez la taille voulue pour N conseillée entre 1024 et 2048\n\n");
		scanf("%d",&Tn_int);
		printf("\n");
		if(Tn_int<=0)
		{
			printf("Erreur de saisie la taille sera 1024\n\n");
			Tn_int = 1024;
		}
		else
		{
			printf("Taille de N = %d\n\n",Tn_int);
		}
		return Tn_int;
}

int main(int argc, char* argv[]){
	unsigned int tab[25] = //tableau contenu les 25 premiers nombres premiers
	{ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};

	int k = sizeof(tab)/sizeof(unsigned int);
	mpz_t n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig;
	mpz_inits(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,NULL);
	int choix;
	mpz_set_ui(e,65537);
	messageATraiter(m);
	printf("Tapez :\n1 pour du chiffrement RSA standard\n2 pour du chiffrement RSA CRT\n3 pour signer le message en standard\n4 pour signer le message en CRT\n\n");
	scanf("%d",&choix);
	printf("\n");
	int taille_bit = taille_N();

	switch(choix)
	{
		case 0: return 0;
		case 1: generation_RSA(n,e,d,taille_bit,tab,k);
				joye_ladder(c,m,e,n);
				joye_ladder(m_obtenu,c,d,n);
				gmp_printf("Cypher_text = %Zx\n\nE : exposant publique = %Zx\n\nN = %Zx\n\n",c,e,n);
				est_egale(m,m_obtenu);
				break;
		case 2: generation_RSA_CRT(n,e,d_p,d_q,I_p,p,q,taille_bit,tab,k);
				joye_ladder(c,m,e,n);
				decrypt_rsa_CRT(m_obtenu,c,d_p,d_q,I_p,p,q);
				gmp_printf("Cypher_text = %Zx\n\nE : exposant publique = %Zx\n\nN = %Zx\n\n",c,e,n);
				est_egale(m,m_obtenu);
				break;
		case 3:	generation_RSA(n,e,d,taille_bit,tab,k);
				joye_ladder(sig,m,d,n);
				joye_ladder(verif_sig,sig,e,n);
				gmp_printf("E : exposant publique pour la vérification = %Zx\n\nN = %Zx\n\nsignature = %Zx\n\n",e,n,verif_sig);
				est_egale(verif_sig,m);
				break;
		case 4: generation_RSA_CRT(n,e,d_p,d_q,I_p,p,q,taille_bit,tab,k);
				sign_CRT(sig,m,p,q,d_p,d_q,I_p,n);
				joye_ladder(verif_sig,sig,e,n);
				gmp_printf("E : exposant publique pour la vérification = %Zx\n\nN = %Zx\n\nsignature = %Zx\n\n",e,n,verif_sig);
				est_egale(verif_sig,m);
				break;
		default : printf("erreur de saisie !\n\n");break;
	}

	mpz_clears(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,NULL);
	return 0;
}

  
