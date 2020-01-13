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

int invert(mpz_t u, mpz_t a, mpz_t n){//version itérative
	mpz_t gcd;
	mpz_init(gcd);
	pgcd_it(gcd,a,n);
	if(mpz_cmp_ui(gcd,1) != 0)
	{
		printf("Elément non inversible !\n");
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
		//gmp_printf("q = %Zd\n",q);
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
		//gmp_printf("u = %Zd\n",u1);
		//gmp_printf("v = %Zd\n",v1);
	}
	mpz_mod(u,u,n);
	mpz_clears(u1,v,v1,r,r1,q,rs,us,vs,tmp1,tmp2,tmp3,NULL);
	return 1;

}

void joye_ladder(mpz_t result, mpz_t m, mpz_t d, mpz_t n){//exponentiation modulaire
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
	//gmp_printf("joye_ladder = %Zd\n",result);
	mpz_clears(tab[0],tab[1],NULL);
}

int test_miller_rabin(mpz_t n, int t){//test de primalité
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

void verif_sign(mpz_t verif_sig, mpz_t m, mpz_t sig, mpz_t e, mpz_t n)//verif_sig = sig^e (mod n) = m
{
	joye_ladder(verif_sig,sig,e,n);
	if(mpz_cmp(verif_sig,m)==0)
	{
		gmp_printf("La signature est valide, on obtient verif_sig = 0x%Zd qui est égal à m\n",verif_sig);
	}
	else
	{
		printf("Problème dans la vérification de signature\n");
		gmp_printf("verif_sig = 0x%Zx",verif_sig);
	}
}
