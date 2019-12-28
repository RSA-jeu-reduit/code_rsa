#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>

//Attention : récursif risqué a cause du nombre d'appel de la fonction (empilement saturé)
//donald knuth : the art of computer programming (analyse du nombre d'appel pour le pgcd en fonction de la taile de a et b)
void pgcd(mpz_t gcd, mpz_t a, mpz_t b, int compteur){ //pour calculer le pgcd 
	//int compteur = 0;
	mpz_t copie_a,copie_b;
	mpz_inits(copie_a,copie_b,NULL);
	mpz_set(copie_a,a);
	mpz_set(copie_b,b);
	if(mpz_cmp_ui(copie_b,0) == 0){
		mpz_set(gcd,copie_a);
	}else{
		mpz_mod(copie_a,copie_a,copie_b);
		compteur++;
		pgcd(gcd,copie_b,copie_a,compteur);
	}
	mpz_clears(copie_a,copie_b,NULL);
	printf("compteur = %d\n",compteur);
}

int AEE(mpz_t u, mpz_t v,mpz_t a, mpz_t n){
	int compteur = 0;
	mpz_t gcd;
	mpz_init(gcd);
	pgcd(gcd,a,n,compteur);
	if (mpz_cmp_ui(gcd,1) != 0){
		printf("Elément non inversible\n");
		return 0;
	}
	mpz_t a1,u1,tmp;
	mpz_inits(a1,u1,tmp,NULL);
	if (mpz_cmp_ui(n,0) == 0){
		mpz_set_ui(u,1);
		mpz_set_ui(v,0);
	}
	else{
		mpz_mod(a1,a,n);//a1 = a mod n
		AEE(u,v,n,a1);//
		mpz_set(u1,u);// u1 <- u
		mpz_set(u,v);//u <- v
		mpz_div(tmp,a,n);
		mpz_mul(tmp,tmp,v);
		mpz_sub(v,u1,tmp);
	}
	mpz_clears(a1,tmp,u1,NULL);
	return 1;
}

void phi_n(mpz_t phi, mpz_t p, mpz_t q){
	mpz_t p_1,q_1;
	mpz_inits(p_1,q_1);
	mpz_sub_ui(p_1,p,1);
	mpz_sub_ui(q_1,q,1);
	mpz_mul(phi,p_1,q_1);	
}
void genere_private_key(mpz_t d, mpz_t e, mpz_t phi){
	mpz_t v;
	mpz_init(v);
	AEE(d,v,e,phi);
	mpz_clear(v);
}

void generation_premier_valable(mpz_t p,mpz_t phip, gmp_randstate_t generateur, mpz_t e, unsigned int long taille_bit){
	int est_premier,premier_entre_eux;
	int compteur = 0;
	mpz_t gcd;
	mpz_init(gcd);
	do{
		mpz_urandomb(p,generateur,taille_bit);
		est_premier = mpz_probab_prime_p(p,4);
		mpz_sub_ui(phip,p,1); 
		pgcd(gcd,e,phip,compteur);
		premier_entre_eux = mpz_cmp_ui(gcd,1);
	}while(premier_entre_eux != 0 || est_premier == 0);
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
	gmp_printf("p = %Zu\nq = %Zu\n\n",p,q);
	gmp_printf("n = 0x%Zx\n\n",n);
	gmp_printf("e = 0x%Zx\n\n",z_e);

	//calcul de d = inv(e) mod phi
	mpz_mul(phi,phip,phiq);
	mpz_invert(d,z_e,phi);
	gmp_printf("d = 0x%Zx\n\n",d);
	mpz_clears(p,q,phip,phiq,phi,NULL);
	gmp_randclear(generateur);
}

/*void chiffrement(mpz_t chiffre, mpz_t message, mpz_t e, mpz_t n){
	puissance(chiffre,message,e,n);
}

void dechiffrement(mpz_t message, mpz_t chiffre, mpz_t d, mpz_t n){
	puissance(message,chiffre,d,n);
}

void dechiffrement_crt(mpz_t chiffre, )*/

void joye_ladder(mpz_t m, mpz_t d, mpz_t n, mpz_t* tab){
	int i,k,bit_d;
	mpz_set_ui(tab[0],1);
	mpz_set(tab[1],m);
	k = mpz_sizeinbase(d,2);
	for(i = 0; i < k; i++){
		bit_d = mpz_tstbit(d,i);
		mpz_mul(tab[1-bit_d],tab[1-bit_d],tab[1-bit_d]);
		mpz_mul(tab[1-bit_d],tab[1-bit_d],tab[bit_d]);
		mpz_mod(tab[1-bit_d],tab[1-bit_d],n);
	}
	gmp_printf("%Zd\n",tab[0]);
}

int main(int argc, char* argv[]){
	mpz_t m,e,d,n,tab[2],gcd;
	mpz_inits(m,e,d,n,tab[0],tab[1],gcd,NULL);
	mpz_set_str(m,"64210519e5344a9c80e70fa7e9a",16);
	printf("ok\n");
	mpz_set_ui(e,3);
	generation_RSA(n,e,d,512);
	joye_ladder(m,d,n,tab);
	gmp_printf("m^d = %Zx\n",tab[0]);
	mpz_powm(m,m,d,n);
	gmp_printf("m^d = %Zx\n",m);
	pgcd(gcd,n,e);

	/*mpz_t a,inv_a,n,v,gcd;
	mpz_inits(a,inv_a,n,v,gcd,NULL);
	mpz_set_ui(n,228);
	mpz_set_ui(v,21);
	mpz_set_ui(inv_a,6);
	mpz_powm(a,inv_a,v,n);
	gmp_printf("%Zd\n",a);*/
	/*mpz_set_ui(a,2345678976578757867);
	mpz_set_ui(n,987674534367567876);
	pgcd(gcd,a,n);
	gmp_printf("%Zd\n",gcd);*/
	//AEE(inv_a,v,a,n);
	//gmp_printf("inv_a = %Zd\n a = %Zd\n",inv_a,a,n);
	//mpz_clears(a,inv_a,n,v,NULL);
	mpz_clears(e,d,n,tab[0],tab[1],m,gcd,NULL);
}
