#include "rsajr.h"


void est_egale(mpz_t m,mpz_t m_obtenu){
	if (mpz_cmp (m_obtenu,m) == 0)
	{
		printf("Déchiffrement ok\n");
	}
	else
	{
		printf("Erreur dans le déchiffrement\n");
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

void crible_opti(mpz_t q, unsigned int b, int k, unsigned long int* tab){
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
	//system("clear");
	mpz_clears(z_residus,NULL);
	for(j = 0; j < k-1; j++){
    	mpz_clear(tab_residus[j]);
    }
    free(tab_residus);
}

void generation_premier_valable(mpz_t p,mpz_t phip, gmp_randstate_t generateur, mpz_t e, unsigned int long taille_bit, int k, unsigned long int* tab){
	mpz_t gcd;
	mpz_init(gcd);

	do{
		mpz_urandomb(p,generateur,taille_bit);

		mpz_setbit(p,taille_bit - 1); //le forcer a être de taille b
		mpz_setbit(p,0);//le forcer à être impaire 

    	crible_opti(p,taille_bit,k,tab);

		mpz_sub_ui(phip,p,1);

		pgcd_it(gcd,e,phip);
	}while(mpz_cmp_ui(gcd,1) != 0);
	mpz_clear(gcd);
}

void generation_RSA(mpz_t n, mpz_t e, mpz_t d, int taille_bit){
	int k = 30;
	mpz_t p,q,phi,phip,phiq;
	mpz_inits(p,q,phi,phip,phiq,NULL);
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	unsigned long int* tab = NULL;
	tab = k_plus_petit_premier(k);
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
	//gmp_printf("p = %Zd\n\n",p);
	//printf("%d\n",mpz_probab_prime_p(p,3));
	//gmp_printf("q = %Zd\n\n",q);
	//printf("%d\n",mpz_probab_prime_p(q,3));
	//calcul de d = inv(e) mod phi
	mpz_mul(phi,phip,phiq);
	mpz_invert(d,e,phi);
	//mpz_mod(d,d,phi);
	//gmp_printf("d = 0x%Zx\n\n",d);
	mpz_inits(p,q,phi,phip,phiq,NULL);
	gmp_randclear(generateur);
}

void generation_RSA_CRT(mpz_t n, mpz_t e, mpz_t d_p, mpz_t d_q, mpz_t I_p, mpz_t p, mpz_t q, unsigned long int taille_bit){
	int k = 30;
	mpz_t phi,phip,phiq;
	mpz_inits(phi,phip,phiq,NULL);
	gmp_randstate_t generateur;
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
	unsigned long int* tab = NULL;
	tab = k_plus_petit_premier(k);
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

void encrypt_rsa(mpz_t c, mpz_t m, mpz_t e, mpz_t n){//c = m^e (mod n)
	joye_ladder(c,m,e,n);
	//mpz_powm(c,m,e,n);
	gmp_printf("cypher = 0x%Zx\n\n",c);
}

void decrypt_rsa(mpz_t m, mpz_t c, mpz_t d, mpz_t n){//m = c^d (mod n)
	joye_ladder(m,c,d,n);
	//mpz_powm(m,c,d,n);
	gmp_printf("m_obtenu = 0x%Zx\n",m);
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
	gmp_printf("m_obtenu = 0x%Zx\n",m);
	mpz_clears(m_p,m_q,NULL);
}

void sign(mpz_t sig, mpz_t m, mpz_t d, mpz_t n)//sig=m^d (mod n)
{
	joye_ladder(sig,m,d,n);
	gmp_printf("signature = 0x%Zx\n", sig);
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


