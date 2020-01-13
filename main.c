#include "rsaJR_crible.h"
#include "rsaJR.h"


int main(int argc, char* argv[]){
	
/*
	int taille_n,taille_phin,e_int;
	printf("Donnez la taille de n voulue, par exemple 1024\n");
	scanf("%d\n",taille_n);

	printf("Donnez la valeur de e voulue, e doit être premier et inférieur à phi(n)= %d, par exemple 65537\n");*/
	//int i;
	long unsigned int taille_bit = 1024; 
	//int k = 30;
	//unsigned long int* tab = NULL;
	//tab = k_plus_petit_premier(k);
	mpz_t n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi;
	mpz_inits(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi,NULL);
	mpz_set_ui(e,65537);
	//gmp_printf("%Zd",e);
	/*FILE* fichier;
	fichier = fopen("timing_projet2_crible.txt","w");
	clock_t start,end;
	double cpu_time_used;
	for (int i = 0; i < 100; i++){
		start = clock();
		generation_RSA(n,e,d,taille_bit);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    	fprintf(fichier,"%d %f\n",i,cpu_time_used);
		printf("%d\n",i);
	}*/
	mpz_set_str(m,"234578435679875467896755435654787687845676743647",10);
	gmp_printf("m = 0x%Zx\n\n",m);
	generation_RSA_CRT(n,e,d_p,d_q,I_p,p,q,1024);
	encrypt_rsa(c,m,e,n);
	decrypt_rsa_CRT(m_obtenu,c,d_p,d_q,I_p,p,q);
	est_egale(m,m_obtenu);
	sign_CRT(sig,m,p,q,d_p,d_q,I_p,n);
	verif_sign(verif_sig,m,sig,e,n);
	//sign(sig,m,d,n);
	/*gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	do{
		generation_premier_valable(p,phip,generateur,e,taille_bit/2,k,tab);
		generation_premier_valable(q,phiq,generateur,e,taille_bit/2 + 1,k,tab);
		mpz_mul(n,p,q);
	}while(mpz_sizeinbase(n,2)!=taille_bit);
	gmp_printf("p = %Zd\n\n",p);
	printf("%d\n",mpz_probab_prime_p(p,3));
	gmp_printf("q = %Zd\n\n",q);
	printf("%d\n",mpz_probab_prime_p(q,3));
	//calcul de d = inv(e) mod phi
	mpz_mul(phi,phip,phiq);
	mpz_invert(d,e,phi);
	//mpz_mod(d,d,phi);
	gmp_printf("d = 0x%Zx\n\n",d);
	//generation_RSA(n,e,d,taille_bit);
	//printf("ok\n");
	mpz_set_str(m,"234578435679875467896755435654787687845676743647",10);
	gmp_printf("m = 0x%Zx\n\n",m);
	encrypt_rsa(c,m,e,n);
	decrypt_rsa(m_obtenu,c,d,n);
	est_egale(m,m_obtenu);*/
	//sign(sig,m,d,n);

	mpz_clears(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi,NULL);


	return 0;
}

/*int main(int argc, char* argv[]){
	int i;
	mpz_t n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig;
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);//initialisation du generateur
	gmp_randseed_ui(generateur,time(NULL));
	mpz_inits(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,NULL);
	mpz_set_ui(e,65537);
	//mpz_set_str(m,"234578435679875467896755435654787687845676743647",10);
	//gmp_printf("m = 0x%Zx\n\n",m);
	//generation_RSA_CRT(n,e,d_p,d_q,I_p,p,q,1024);
	//encrypt_rsa(c,m,e,n);
	//decrypt_rsa_CRT(m_obtenu,c,d_p,d_q,I_p,p,q);
	//est_egale(m,m_obtenu);
	//sign(sig,m,d,n);
	FILE* fichier;
	fichier = fopen("timing_projet2.txt","w");
	clock_t start,end;
	double cpu_time_used;
	for (i = 0; i < 100; i++){
		start = clock();
		generation_RSA(n,e,d,1024);
		end = clock();
		cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    	fprintf(fichier,"%d %f\n",i,cpu_time_used);
		printf("%d\n",i);
	}
	/*generation_RSA(n,e,d,1024);
	mpz_set_str(m,"234578435679875467896755435654787687845676743647",10);
	gmp_printf("m = 0x%Zx\n\n",m);
	encrypt_rsa(c,m,e,n);
	decrypt_rsa(m_obtenu,c,d,n);
	est_egale(m,m_obtenu);*/
	//sign(sig,m,d,n);
	/*gmp_randclear(generateur);
	mpz_clears(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,NULL);


	return 0;
}*/
