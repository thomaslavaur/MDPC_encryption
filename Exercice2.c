#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define w 39//39
#define n 4813//4813
#define e_weight 2*w
#define T 26//26
#define message "Lorem ipsum dolor sit amet,"


struct polynome
{
	int dimension;
	int* valeurs;
};

typedef struct polynome polynome;


struct cle_privee
{
	polynome h0;
	polynome h1;
};

typedef struct cle_privee cle_privee;


struct cle
{
	cle_privee sk;
	polynome pk;
};

typedef struct cle cle;


struct chiffre
{
	polynome c0;
	polynome c1;
};

typedef struct chiffre chiffre;

struct error
{
	polynome e0;
	polynome e1;
};

typedef struct error error;


int val_pol(polynome P, int indice)
{
	return(*(P.valeurs+indice));
};


int* pt(polynome P, int indice)
{
	return(P.valeurs+indice);
};


int val_mat(polynome P, int ligne, int colonne)
{
	return(*(P.valeurs+((colonne - ligne + P.dimension) % P.dimension)));
};


void destruction_polynome(polynome P)
{
	free(P.valeurs);
	P.valeurs = NULL;
};


void destruction_cle(cle K)
{
	destruction_polynome(K.sk.h0);
	destruction_polynome(K.sk.h1);
	destruction_polynome(K.pk);
};


int poids_hamming(polynome P)
{
	int poids = 0;
	for(int i = 0; i < P.dimension; i++)
	{
		if(val_pol(P,i) != 0)
		{
			poids++;
		}
	}
	return poids;
};


polynome generation_polynome(int dimension,int poids)
{
	polynome P;
	P.dimension = dimension;
	P.valeurs = (int*) malloc(P.dimension*sizeof(int));

	int j;
	int memoire;
	for(int i = 0; i < dimension; i++)
	{
		if(i < poids)
		{
			*pt(P,i) = 1;
		}
		else
		{
			*pt(P,i) = 0;
		}
	}

	for(int i = dimension-1; i > 0; i--)
	{
		j = rand() % (i+1);
		memoire = val_pol(P,i);
		*pt(P,i) = val_pol(P,j);
		*pt(P,j) = memoire;
	}

	return(P);
};


polynome multiplication_polynome(polynome A, polynome B)
{
	polynome resultat;
	resultat.dimension = A.dimension;
	resultat.valeurs = (int*) malloc(resultat.dimension*sizeof(int));

	for(int i = 0; i < resultat.dimension; i++)
	{
		*pt(resultat,i) = 0;
		for(int j = 0; j < resultat.dimension; j++)
		{
			*pt(resultat,i) = val_pol(resultat,i) ^ val_pol(A,j)*val_pol(B,(i-j + B.dimension) % B.dimension); 
		}
	}

	return resultat;
};


polynome inversion_polynome(polynome P)
{
	int* copie;
	int* inverse;
	polynome resultat;

	copie = (int*) malloc(P.dimension*P.dimension*sizeof(int));
	inverse = (int*) malloc(P.dimension*P.dimension*sizeof(int));
	resultat.dimension = P.dimension;
	resultat.valeurs = (int*) malloc(P.dimension*sizeof(int));

	for(int i = 0; i < P.dimension; i++)
	{
		for(int j = 0; j < P.dimension; j++)
		{
			if(i == j)
			{
				*(inverse + i*P.dimension + j) = 1;
			}
			else
			{
				*(inverse + i*P.dimension + j) = 0;
			}
			*(copie + i*P.dimension + j) = val_mat(P,i,j);
		}
	}

	int l;
	int memoire;
	int i = 0;
	int j = 0;

	while((i < P.dimension) &&(j < P.dimension))
	{
		l = i;
		while(*(copie + l*P.dimension + j) != 1)
		{
			l++;
			if(l == P.dimension)
			{
				free(copie);
				free(inverse);
				destruction_polynome(resultat);
				resultat.dimension = 0;
				return(resultat);
			}
		}

		for(int h = 0; h < P.dimension; h++)
		{
			memoire = *(copie + i*P.dimension + h);
			*(copie + i*P.dimension + h) = *(copie + l*P.dimension + h);
			*(copie + l*P.dimension + h) = memoire;
			memoire = *(inverse + i*P.dimension + h);
			*(inverse + i*P.dimension + h) = *(inverse + l*P.dimension + h);
			*(inverse + l*P.dimension + h) = memoire;
		}

		for(int ligne = 0; ligne < P.dimension; ligne++)
		{
			if((ligne != i) && (*(copie + ligne*P.dimension + j) == 1))
			{
				for(int colonne = 0; colonne < P.dimension; colonne++)
				{
					*(copie + ligne*P.dimension + colonne) = (*(copie + ligne*P.dimension + colonne) ^ *(copie + i*P.dimension + colonne));
					*(inverse + ligne*P.dimension + colonne) = (*(inverse + ligne*P.dimension + colonne) ^ *(inverse + i*P.dimension + colonne));
				}
			}
		}

		i++;
		j++;
	}

	for(int i = 0; i < P.dimension; i++)
	{
		*pt(resultat,i) = *(inverse + i);
	}

	free(copie);
	free(inverse);
	return(resultat);
};


cle keygen(int longueur, int poids)
{
	cle resultat;
	
	resultat.sk.h1 = generation_polynome(longueur,poids);

	polynome inverse;

	do
	{
		resultat.sk.h0 = generation_polynome(longueur,poids);

		inverse = inversion_polynome(resultat.sk.h0);

		if(inverse.dimension == 0)
		{
			destruction_polynome(resultat.sk.h0);
		}
	}while(inverse.dimension == 0);

	resultat.pk = multiplication_polynome(resultat.sk.h1,inverse);

	destruction_polynome(inverse);

	return(resultat);
};


polynome hash_polynomes(polynome e0, polynome e1)
{
	char* mot;
	mot = (char*) malloc((e0.dimension+e1.dimension)*sizeof(char));

	for(int i = 0; i < e0.dimension; i++)
	{
		sprintf(mot+i, "%d", val_pol(e0,i));
	}
	for(int i = 0; i < e1.dimension; i++)
	{
		sprintf(mot+e0.dimension+i, "%d", val_pol(e1,i));
	}

	char commande[2*n+29] = "echo -n ";
	strcat(commande,mot);
	strcat(commande," | openssl sha256 -r");

	char byte[2];
	FILE* hash = popen(commande,"r");

	polynome resultat;
	resultat.dimension = 32;
	resultat.valeurs = (int*) malloc(resultat.dimension*sizeof(int));

	for(int i = 0; i < 32; i++)
	{	
		fscanf(hash," %2s",byte);
		*pt(resultat,i) = strtol(byte,NULL,16);
	}

	pclose(hash);								
	free(mot);

	return resultat;
};


polynome XOR(polynome P, polynome hash)
{
	polynome resultat;
	resultat.dimension = P.dimension;
	resultat.valeurs = (int*) malloc(resultat.dimension*sizeof(int));

	for(int i = 0; i < P.dimension; i++)
	{
		*pt(resultat,i) = val_pol(P,i) ^ val_pol(hash,i);
	}

	return resultat;
};


chiffre Chiffrement(polynome h)
{
	chiffre c;

	if(strlen(message) > 32)
	{
		c.c0.dimension = 0;
		c.c1.dimension = 0;
		return c;
	}

	int poidse0 = rand() % e_weight;
	polynome e0 = generation_polynome(n,poidse0);
	polynome e1 = generation_polynome(n,e_weight-poidse0);

	c.c1 = multiplication_polynome(h,e1);
	for(int i = 0; i < c.c1.dimension; i++)
	{
		*pt(c.c1,i) = val_pol(c.c1,i) ^ val_pol(e0,i);
	}

	polynome hash = hash_polynomes(e0,e1);

	destruction_polynome(e0);
	destruction_polynome(e1);


	c.c0.dimension = strlen(message);
	c.c0.valeurs = (int*) malloc(c.c0.dimension*sizeof(int));
	for(int i = 0; i < c.c0.dimension; i++)
	{
		*pt(c.c0,i) = (int) message[i] ^ val_pol(hash,i);
	}
	destruction_polynome(hash);

	return c;
};


error bitflip(polynome h0, polynome h1, polynome s, int seuil, int t)
{
	polynome u, v;
	u.dimension = n;
	v.dimension = n;
	u.valeurs = (int*) malloc(n*sizeof(int));
	v.valeurs = (int*) malloc(n*sizeof(int));

	for(int i = 0; i < n; i++)										
	{																	
		*pt(u,i) = 0;												
		*pt(v,i) = 0;															
	}		


	polynome syndrome;
	syndrome.dimension = n;
	syndrome.valeurs = (int*) malloc(n*sizeof(int));

	for(int i = 0; i < n; i++)
	{
		*pt(syndrome,i) = val_pol(s,i);
	}


	polynome sum;
	sum.dimension = 2*n;
	sum.valeurs = (int*) malloc(2*n*sizeof(int));

	polynome flipped_pos;
	flipped_pos.dimension = 2*n;
	flipped_pos.valeurs = (int*) malloc(2*n*sizeof(int));


	while(((poids_hamming(u) != t) || (poids_hamming(v) != t)) && (poids_hamming(syndrome) != 0))
	{
		for(int i = 0; i < 2*n; i++)
		{
			*pt(sum,i) = 0;
			for(int j = 0; j < n; j++)
			{
				if(i < n)
				{
					*pt(sum,i) = val_pol(sum,i) + (val_pol(syndrome,j) * val_mat(h0,i,j));
				}
				else
				{
					*pt(sum,i) = val_pol(sum,i) + (val_pol(syndrome,j) * val_mat(h1,i-n,j));
				}
			}
		}

		for(int i = 0; i < 2*n; i++)
		{
			if(val_pol(sum,i) >= seuil)
			{
				*pt(flipped_pos,i) = 1;
			}
			else
			{
				*pt(flipped_pos,i) = 0;
			}
		}


		for(int i = 0; i < n; i++)
		{
			*pt(u,i) = val_pol(u,i) ^ val_pol(flipped_pos,i);
			*pt(v,i) = val_pol(v,i) ^ val_pol(flipped_pos, n+i);
		}


		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < 2*n; j++)
			{
				if(j < n)
				{
					*pt(syndrome,i) = val_pol(syndrome,i) ^ (val_pol(flipped_pos,j) * val_mat(h0,j,i));
				}
				else
				{
					*pt(syndrome,i) = val_pol(syndrome,i) ^ (val_pol(flipped_pos,j) * val_mat(h1,j-n,i));
				}			
			}
		}
	}



	for(int i = 0; i < n; i++)
	{
		*pt(syndrome,i) = val_pol(s,i);
		for(int j = 0; j < 2*n; j++)
		{
			if(j < n)
			{
				*pt(syndrome,i) = val_pol(syndrome,i) ^ (val_pol(u,j) * val_mat(h0,j,i));
			}
			else
			{
				*pt(syndrome,i) = val_pol(syndrome,i) ^ (val_pol(v,j-n) * val_mat(h1,j-n,i));
			}
		}
	}


	destruction_polynome(sum);
	destruction_polynome(flipped_pos);
	error resultat;

	if(poids_hamming(syndrome) != 0)
	{
		destruction_polynome(syndrome);
		destruction_polynome(u);
		destruction_polynome(v);
		resultat.e0.dimension = 0;
		resultat.e1.dimension = 0;
		return resultat;
	}
	else
	{
		destruction_polynome(syndrome);
		resultat.e0.dimension = u.dimension;
		resultat.e1.dimension = v.dimension;
		resultat.e0.valeurs = u.valeurs;
		resultat.e1.valeurs = v.valeurs;
		return(resultat);
	}
};


void Dechiffrement(chiffre c, cle_privee sk)
{

	polynome s = multiplication_polynome(sk.h0,c.c1);

	error e;

	e = bitflip(sk.h0,sk.h1,s,T,w);

	if(e.e0.dimension == 0)
	{
		printf("\nERREUR BITFLIP\n");
		return EXIT_FAILURE;
	}

	char* dechiffre;
	dechiffre = (char*) malloc((c.c0.dimension+1)*sizeof(char));

	polynome hash = hash_polynomes(e.e0,e.e1);

	polynome dechiffrage = XOR(c.c0,hash);

	for(int i = 0; i < dechiffrage.dimension; i++)
	{
		*(dechiffre + i) = (char) val_pol(dechiffrage,i);
	}

	*(dechiffre + c.c0.dimension) = '\0';



	printf("\nLe déchiffrer est : %s\n\n",dechiffre);

	destruction_polynome(s);
	destruction_polynome(e.e0);
	destruction_polynome(e.e1);
	free(dechiffre);
	destruction_polynome(hash);
	destruction_polynome(dechiffrage);
};


int main()
{
	srand(time(NULL));
	float temps1 = 0;
	float temps2 = 0;
	float temps3 = 0;
	float temps4 = 0;

	cle MDPC;

	temps1 = clock();
	MDPC = keygen(n,w);
	temps2 = clock();

	chiffre c = Chiffrement(MDPC.pk);
	temps3 = clock();

	Dechiffrement(c,MDPC.sk);
	temps4 = clock();

	printf("\nLe message à crypter est : %s",message);

	temps1 = (float)(temps2-temps1)/CLOCKS_PER_SEC;
	temps2 = (float)(temps3-temps2)/CLOCKS_PER_SEC;
	temps3 = (float)(temps4-temps3)/CLOCKS_PER_SEC;

	printf("\n\nTemps de la génération des clés : %f s\nTemps du chiffrement d'un bloc : %f s\nTemps du déchiffrement d'un bloc : %f s\nTemps total : %f s\n\n",temps1,temps2,temps3,temps1+temps2+temps3);

	destruction_polynome(c.c0);
	destruction_polynome(c.c1);
	destruction_cle(MDPC);
}