#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#define n 1000
#define k 500
#define t 10


struct matrice
{
	unsigned int colonnes;
	unsigned int lignes;
	unsigned int *valeurs;
};

typedef struct matrice matrice;


unsigned int val(matrice A, unsigned int ligne, unsigned int colonne)
{
	return(*(A.valeurs+(ligne*A.colonnes)+colonne));
};

unsigned int* pt(matrice A, unsigned int ligne, unsigned int colonne)
{
	return(A.valeurs+(ligne*A.colonnes)+colonne);
};

unsigned int bit_generator()
{
	static unsigned int compteur = 15;
	static unsigned int random;
	unsigned int resultat;
	if(compteur == 15)
	{
		compteur = 0;
		random = rand();
	}
	resultat = 1 & (random >> compteur);
	compteur ++;
	return(resultat);
};


void destruction_matrice(matrice M)
{
	free(M.valeurs);
	M.valeurs = NULL;
};


unsigned int poids_hamming(matrice mot)
{
	unsigned int poids = 0;
	for(unsigned int i = 0; i < mot.colonnes; i++)
	{
		if(val(mot,0,i) == 1)
		{
			poids++;
		}
	}
	return poids;
};


matrice produit_matrice(matrice M, matrice N)
{
	if(M.colonnes != N.lignes)
	{
		printf("Erreur de produit : M n'a pas le même nombre de colonne que N a de lignes");
		return(M);
	}
	else
	{
		matrice resultat;
		resultat.valeurs = (unsigned int*) malloc(M.lignes*N.colonnes*sizeof(unsigned int));
		resultat.lignes = M.lignes;
		resultat.colonnes = N.colonnes;
		for(unsigned int i = 0; i< M.lignes ; i++)
		{
			for(unsigned int j = 0 ; j < N.colonnes ; j++)
			{
				*pt(resultat,i,j) = 0;
				for(unsigned int l = 0; l < M.colonnes ; l++)
				{
					*pt(resultat,i,j) = val(resultat,i,j) ^ val(M,i,l)*val(N,l,j);
				}
			}
		}
		return(resultat);
	}
};


matrice inverse_matrice(matrice A)
{
	matrice M;
	M.lignes = A.lignes;
	M.colonnes = A.colonnes;
	M.valeurs = (unsigned int*) malloc(M.lignes*M.colonnes*sizeof(unsigned int));

	matrice inverse;
	inverse.lignes = A.lignes;
	inverse.colonnes = A.colonnes;
	inverse.valeurs = (unsigned int*) malloc(inverse.lignes*inverse.colonnes*sizeof(unsigned int));

	for(unsigned int i = 0; i < A.lignes; i++)
	{
		for(unsigned int j = 0; j < A.colonnes; j++)
		{
			if(i == j)
			{
				*pt(inverse,i,j) = 1;
			}
			else
			{
				*pt(inverse,i,j) = 0;
			}
			*pt(M,i,j) = val(A,i,j);
		}
	}

	unsigned int l;
	unsigned int memoire;
	for(unsigned int colonne = 0; colonne < A.lignes; colonne ++)
	{
		l = colonne;
		while(val(M,l,colonne) != 1)
		{
			l++;
		}
		if(l != colonne)
		{
			for(unsigned int temp = 0; temp < A.lignes; temp++)
			{
				memoire = val(M,l,temp);
				*pt(M,l,temp) = val(M,colonne,temp);
				*pt(M,colonne,temp) = memoire;
				memoire = val(inverse,l,temp);
				*pt(inverse,l,temp) = val(inverse,colonne,temp);
				*pt(inverse,colonne,temp) = memoire;
			}
		}
		for(unsigned int ligne = 0; ligne < A.lignes; ligne++)
		{
			if((ligne != colonne) && (val(M,ligne,colonne) == 1))
			{
				for(unsigned int temp = 0; temp < A.lignes; temp++)
				{
					*pt(M,ligne,temp) = val(M,ligne,temp) ^ val(M,colonne,temp);
					*pt(inverse,ligne,temp) = val(inverse,ligne,temp) ^ val(inverse,colonne,temp);
				}
			}
		}
	}
	destruction_matrice(M);
	return(inverse);
};


matrice transposition_matrice(matrice A)
{
	matrice T;
	T.lignes = A.colonnes;
	T.colonnes = A.lignes;
	T.valeurs = (unsigned int*) malloc(T.lignes*T.colonnes*sizeof(unsigned int));

	for(unsigned int i = 0; i < T.lignes; i++)
	{
		for(unsigned int j = 0; j < T.colonnes; j++)
		{
			*pt(T,i,j) = val(A,j,i);
		}
	}
	return(T);
};


matrice generation_permutation(unsigned int dimension)
{
	matrice P;
	P.lignes = dimension;
	P.colonnes = dimension;
	P.valeurs = (unsigned int*) malloc(P.lignes*P.colonnes*sizeof(unsigned int));

	matrice liste;
	liste.lignes = 1;
	liste.colonnes = dimension;
	liste.valeurs = (unsigned int*) malloc(liste.lignes*liste.colonnes*sizeof(unsigned int));
	unsigned int j;
	unsigned int memoire;
	for(unsigned int i = 0; i < dimension; i++)
	{
		*pt(liste,0,i) = i;
	}
	for(unsigned int i = dimension-1; i > 0; i--)
	{
		j = rand() % (i+1);
		memoire = val(liste,0,i);
		*pt(liste,0,i) = val(liste,0,j);
		*pt(liste,0,j) = memoire;
	}

	for(unsigned int i = 0; i < dimension; i++)
	{
		for(unsigned int j = 0; j < dimension; j++)
		{
			if(j == val(liste,0,i))
			{
				*pt(P,i,j) = 1;
			}
			else
			{
				*pt(P,i,j) = 0;
			}
		}
	}

	destruction_matrice(liste);
	return(P);
};


matrice Gauss_jordan(matrice A) // renvoie U tel que UA = (I|H)
{
	matrice U;
	U.lignes = A.lignes;
	U.colonnes = A.lignes;
	U.valeurs = (unsigned int*) malloc(U.lignes*U.colonnes*sizeof(unsigned int));

	matrice M;
	M.lignes = A.lignes;
	M.colonnes = A.colonnes;
	M.valeurs = (unsigned int*) malloc(M.lignes*M.colonnes*sizeof(unsigned int));

	for(unsigned int i = 0; i < U.lignes; i++)
	{
		for(unsigned int j = 0; j < U.colonnes; j++)
		{
			if(i == j)
			{
				*pt(U,i,j) = 1;
			}
			else
			{
				*pt(U,i,j) = 0;
			}
		}
	}

	for(unsigned int i = 0; i < A.lignes; i++)
	{
		for(unsigned int j = 0; j < A.colonnes; j++)
		{
			*pt(M,i,j) = val(A,i,j);
		}
	}
	
	unsigned int l;
	unsigned int memoire;
	unsigned int i = 0;
	unsigned int j = 0;

	while((i < M.lignes) &&(j < M.colonnes))
	{
		l = i;
		while(val(M,l,j) != 1)
		{
			l++;
			if(l == M.lignes)
			{
				destruction_matrice(U);
				destruction_matrice(M);
				U.lignes = 0;
				U.colonnes = 0;
				return(U);
			}
		}

		for(unsigned int h = 0; h < M.colonnes; h++)
		{
			memoire = val(M,i,h);
			*pt(M,i,h) = val(M,l,h);
			*pt(M,l,h) = memoire;
		}
		for(unsigned int h = 0; h < U.colonnes; h++)
		{
			memoire = val(U,i,h);
			*pt(U,i,h) = val(U,l,h);
			*pt(U,l,h) = memoire;
		}

		for(unsigned int ligne = 0; ligne < M.lignes; ligne++)
		{
			if((ligne != i) && (val(M,ligne,j) == 1))
			{
				for(unsigned int colonne = 0; colonne < M.colonnes; colonne++)
				{
					*pt(M,ligne,colonne) = val(M,ligne,colonne) ^ val(M,i,colonne);
				}
				for(unsigned int colonne = 0; colonne < U.colonnes; colonne++)
				{
					*pt(U,ligne,colonne) = val(U,ligne,colonne) ^ val(U,i,colonne);
				}
			}
		}

		i++;
		j++;
	}
	destruction_matrice(M);
	return(U);
};


matrice Prange_ISD(matrice H, matrice S, unsigned int w)
{
	unsigned int poids;
	matrice P;
	matrice HP;
	matrice U;
	matrice Ut;
	matrice SUt;
	matrice resultat;
	matrice SUt0;
	matrice P_inverse;

	do
	{
		do
		{
			P = generation_permutation(H.colonnes);
			HP = produit_matrice(H,P);
			U = Gauss_jordan(HP);
			if((U.lignes == 0) && (U.colonnes == 0))
			{
				destruction_matrice(P);
				destruction_matrice(HP);
			}
		}while((U.lignes == 0) && (U.colonnes == 0));
		Ut = transposition_matrice(U);
		SUt = produit_matrice(S,Ut);
		poids = poids_hamming(SUt);
		if( poids < w)
		{
			destruction_matrice(P);
			destruction_matrice(HP);
			destruction_matrice(U);
			destruction_matrice(Ut);
			destruction_matrice(SUt);
		}
	}while(poids < w);

	SUt0.lignes = 1;
	SUt0.colonnes = H.colonnes;
	SUt0.valeurs = (unsigned int*) malloc(SUt0.lignes*SUt0.colonnes*sizeof(unsigned int));

	for(unsigned int i = 0; i < SUt.colonnes; i++)
	{
		*pt(SUt0,0,i) = val(SUt,0,i);
	}
	for(unsigned int i = SUt.colonnes; i < SUt0.colonnes; i++)
	{
		*pt(SUt0,0,i) = 0;
	}

	P_inverse = inverse_matrice(P);
	resultat = produit_matrice(SUt0,P_inverse);

	destruction_matrice(P);
	destruction_matrice(HP);
	destruction_matrice(U);
	destruction_matrice(Ut);
	destruction_matrice(SUt);
	destruction_matrice(SUt0);
	destruction_matrice(P_inverse);

	return(resultat);
};


matrice matrice_aleatoire(unsigned int ligne, unsigned int colonne)
{
	matrice resultat;
	resultat.lignes = ligne;
	resultat.colonnes = colonne;
	resultat.valeurs = (unsigned int*) malloc(resultat.lignes*resultat.colonnes*sizeof(unsigned int));

	for(unsigned int i = 0; i < ligne; i++)
	{
		for(unsigned int j = 0; j < colonne; j++)
		{
			*pt(resultat,i,j) = bit_generator();
		}
	}

	return(resultat);
};


int main()
{
	srand(time(NULL));
	float temps1 = 0;
	float temps2 = 0;

	printf("\n%d",i);
	matrice H = matrice_aleatoire(n-k,n);

	matrice s = matrice_aleatoire(1,n-k);

	unsigned int w = t;

	temps1 = clock();
	matrice resultat = Prange_ISD(H,s,w);
	temps2 = clock();
	total1 = (float)(temps2-temps1)/CLOCKS_PER_SEC;

	destruction_matrice(resultat);
	destruction_matrice(s);
	destruction_matrice(H);

	printf("\n\nTemps d'exection : %f\n\n",temps1);
	return(0);
}