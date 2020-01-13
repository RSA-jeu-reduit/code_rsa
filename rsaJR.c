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

void generation_premier_valable(mpz_t p,mpz_t phip, gmp_randstate_t generateur, mpz_t e, unsigned int long taille_bit){
	mpz_t gcd;
	mpz_init(gcd);
	//printf("io");
	do{
		mpz_urandomb(p,generateur,taille_bit);
		mpz_sub_ui(phip,p,1); 
		pgcd_it(gcd,e,phip);
	}while(mpz_cmp_ui(gcd,1) != 0 || test_miller_rabin(p,3) == 0);
	mpz_clear(gcd);
}

void generation_RSA(mpz_t n, mpz_t z_e,mpz_t d,unsigned long int k){
	mpz_t p,q,phip,phiq,phi;
	gmp_randstate_t generateur;
	unsigned long int taille_bit;
	mpz_inits(p,q,phip,phiq,phi,NULL);
	gmp_randinit_default(generateur);
	gmp_randseed_ui(generateur,time(NULL));
		do{ 
		if (k%2 == 0){//distinction des cas où k est pair ou impair
			generation_premier_valable(p,phip,generateur,z_e,k/2);
			generation_premier_valable(q,phiq,generateur,z_e,k/2);
		}
		else{
			generation_premier_valable(p,phip,generateur,z_e,k/2 + 1);
			generation_premier_valable(q,phiq,generateur,z_e,k/2);
		}
		mpz_mul(n,p,q);
		taille_bit = mpz_sizeinbase(n,2);
	}while (taille_bit != k);
	//gmp_printf("phip = %Zu\nphiq = %Zu\n\n",phip,phiq);
	//gmp_printf("p = %Zu\nq = %Zu\n\n",p,q);
	//gmp_printf("n = 0x%Zx\n\n",n);
	//gmp_printf("e = 0x%Zx\n\n",z_e);

	//calcul de d = inv(e) mod phi
	mpz_mul(phi,phip,phiq);
	invert(d,z_e,phi);
	//mpz_mod(d,d,phi);
	//gmp_printf("d = 0x%Zx\n\n",d);
	mpz_clears(p,q,phip,phiq,phi,NULL);
	gmp_randclear(generateur);
}

void generation_RSA_CRT(mpz_t n, mpz_t z_e, mpz_t d_p, mpz_t d_q, mpz_t I_p, mpz_t p, mpz_t q, unsigned long int k){
	mpz_t phip,phiq,phi;
	gmp_randstate_t generateur;
	unsigned long int taille_bit;
	mpz_inits(phip,phiq,phi,NULL);
	gmp_randinit_default(generateur);//initialisation du generateur
	gmp_randseed_ui(generateur,time(NULL));
		do{
		if (k%2 == 0){
			generation_premier_valable(p,phip,generateur,z_e,k/2);
			generation_premier_valable(q,phiq,generateur,z_e,k/2);
		}
		else{
			generation_premier_valable(p,phip,generateur,z_e,k/2 + 1);
			generation_premier_valable(q,phiq,generateur,z_e,k/2);
		}
		mpz_mul(n,p,q);
		taille_bit = mpz_sizeinbase(n,2);
	}while (taille_bit != k);
	invert(d_p,z_e,phip);
	invert(d_q,z_e,phiq);
	invert(I_p,p,q);
	/*gmp_printf("p = %Zu\nq = %Zu\n\n",p,q);
	gmp_printf("n = 0x%Zx\n\n",n);
	gmp_printf("d_p = 0x%Zx\nd_q = %Zu\n\n",d_p,d_q);
	gmp_printf("I_p = 0x%Zx\n\n",I_p);*/
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

void verif_sign(mpz_t verif_sig, mpz_t m, mpz_t sig, mpz_t e, mpz_t n)//verif_sig = sig^e (mod n) = m
{
	joye_ladder(verif_sig,sig,e,n);
	if(mpz_cmp(verif_sig,m)==0)
	{
		gmp_printf("La signature est valide, on obtient verif_sig = 0x%Zx qui est égal à m\n",verif_sig);
	}
	else
	{
		printf("Problème dans la vérification de signature\n");
		gmp_printf("verif_sig = 0x%Zx",verif_sig);
	}
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

