int main(int argc, char* argv[]){
	mpz_t n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi;
	mpz_inits(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi,NULL);
	int choix;

	messageATraiter(m);
	printf("Tapez :\n1 pour du chiffrement RSA standard\n2 pour du chiffrement RSA CRT\n3 pour signer le message en standard\n4 en CRT\n0 pour arrêter\n");
	scanf("%d\n",&choix);
	int taille_bit = taille_N();
	choixE(e,taille_bit);

	switch(choix)
	{
		case 0: return 0;
		case 1: generation_RSA(n,e,d,taille_bit);
				joye_ladder(c,m,e,n);
				joye_ladder(m_obtenu,c,d,n);
				est_egale(m,m_obtenu);break;
		case 2: generation_RSA_CRT(n,e,d_p,d_q,I_p,p,q,taille_bit);
				joye_ladder(c,m,e,n);
				decrypt_rsa_CRT(m_obtenu,c,d_p,d_q,I_p,p,q);
				est_egale(m,m_obtenu);break;
		case 3:	joye_ladder(sig,m,d,n);
				joye_ladder(verif_sig,sig,e,n);break;
				est_egale(verif_sig,m);break;
		case 4: sign_CRT(sig,m,p,q,d_p,d_q,I_p,n);
				joye_ladder(verif_sig,sig,e,n);
				est_egale(verif_sig,m);break;
		default : printf("erreur de saisie\n");break;
	}

	mpz_clears(n,e,d,m,m_obtenu,c,p,q,d_p,d_q,I_p,sig,verif_sig,phip,phiq,phi,NULL);
	return 0;
}

   
