#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define P 4
#define PP 2 //liczba procesorów PP*PP - algorytm dla P^2 procesorów
#define N 2016 //rozmiar tablic

float A_Odczyt[N][N], B_Odczyt[N][N], C_Row[N][N], C_Sek[N][N];
float A[N / PP][N / PP], B[N / PP][N / PP], C[N / PP][N / PP];
float AA[N / PP][N / PP], BB[N / PP][N / PP];
float A_Temp[N / PP][N / PP], B_Temp[N / PP][N / PP], C_Temp[N / PP][N / PP];

float(*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];

double startwtime1, startwtime2, startwtime3, endwtime;

void print_matrix(float matrix[N][N])
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++) 
        {
            printf("%.1f ", matrix[i][j]);
        }
        printf("\n");
    }
}

int przesun_prawo(int id_procesu, int ile_przesuniec)
{
    if (ile_przesuniec == 0)
    {
        return id_procesu;
    }

    int wynik = id_procesu;

    for (int i = 0; i < ile_przesuniec; i++)
    {
        wynik = ((wynik + 1) % PP) + ((wynik / PP) * PP);
    }

    return wynik;
}

int przesun_dol(int id_procesu, int ile_przesuniec)
{
    if (ile_przesuniec == 0)
    {
        return id_procesu;
    }

    int wynik = id_procesu;

    for (int i = 0; i < ile_przesuniec; i++)
    {
        wynik = (wynik + PP) % P;
    }

    return wynik;
}

int modulo(int x)
{
    int reszta = x % PP;
    if (reszta < 0)
    {
        return reszta + PP;
    }
    else
    {
        return reszta;
    }
}

int main(int argc, char** argv)
{

	FILE* plik;
	FILE* plik_out;

	int wiersz, kolumna, mod = 0;
	int data_received = -1;
	int tag = 101;
	int koniec;

    int ilosc_procesow;
    int id_procesu;


	MPI_Status statSend[2], statRecv[2], statRecvZbierz[P];
	MPI_Request reqSend[2], reqRecv[2], reqSendZbierz[P], reqRecvZbierz[P];

	MPI_Init(0, 0);
	MPI_Comm_rank(MPI_COMM_WORLD, &id_procesu);
	MPI_Comm_size(MPI_COMM_WORLD, &ilosc_procesow);

	int Lewo = modulo(id_procesu - 1) + ((id_procesu / PP) * PP);
	int Prawo = ((id_procesu + 1) % PP) + ((id_procesu / PP) * PP);
	int Dol = (id_procesu + PP) % P;
    int Gora = (id_procesu + P - PP) % P;

	if (id_procesu == 0)
    {
        printf("--------------------------------------------------------------------------------\n");
		printf("Metoda Cannona dla tablicy %d x %d elementow \n", N, N);
    }

	if (id_procesu == 0) 
    {
        startwtime1 = MPI_Wtime();//czas w sekundach
    }

	//wczytanie danych przez proces rank=0
	if (id_procesu == 0)
	{
		plik = fopen("liczby.txt", "r");
		if (plik == NULL)
		{
			printf("Blad otwarcia pliku \"liczby.txt\"\n");
			koniec = 1;
			MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Finalize();
			exit(0);
		}
		else 
        {
			koniec = 0;
			MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if (koniec) 
        { 
            MPI_Finalize(); 
            exit(0); 
        }
	}

	if (ilosc_procesow != P) 
    {
		if (id_procesu == 0) 
        {   
            printf("wywolano obliczenia iloczynu macierzy metoda cannona na %d procesach - uruchom mpiexec -n %d matrixmult\n", ilosc_procesow, P);
        }

        MPI_Finalize();
        exit(0);
	}

	if (id_procesu == 0)
	{
        for (int i = 0; i < N ; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (fscanf(plik, "%f", &A_Odczyt[i][j]) != 1)
                {
                    printf("Blad podczas wczytywania A.\n");
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
                else
                {
                    //printf("Wczytano A[%d][%d] = %f\n", i, j, A_Odczyt[i][j]);
                }
            }
        }

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (fscanf(plik, "%f", &B_Odczyt[i][j]) != 1)
                {
                    printf("Blad podczas wczytywania B.\n");
                    MPI_Finalize();
                    exit(EXIT_FAILURE);
                }
                else
                {
                    //printf("Wczytano B[%d][%d] = %f\n", i, j, B_Odczyt[i][j]);
                }
            }
        }


        for(int i = 0; i < N / PP; i++)
        {
  			for(int j = 0 ; j < N / PP; j++)
            {
  				A[i][j] = A_Odczyt[i][j];
  				B[i][j] = B_Odczyt[i][j];
  			}
  		}

        for(int i = 1; i < P; i++)
        {
  			for(int j = 0; j < N/PP; j++)
            {
  				for(int k = 0; k < N/PP; k++)
                {
					A_Temp[j][k] = A_Odczyt[(przesun_prawo(i, i / PP) / PP) * (N / PP) + j][(przesun_prawo(i, i / PP) % PP) * (N / PP) + k];
					
					B_Temp[j][k] = B_Odczyt[(przesun_dol(i, i % PP) / PP) * (N / PP) + j][(przesun_dol(i, i % PP) % PP) * (N / PP) + k];
				}
			}
			
			MPI_Isend(A_Temp, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[0]);
			MPI_Isend(B_Temp, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[1]);
		}

	}
	else
	{
		MPI_Irecv(A, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[0]);
		MPI_Irecv(B, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[1]);
		
		MPI_Waitall(2, reqRecv, statRecv);
	}

	for (int i = 0; i < N / PP; i++)
    {
		for (int j = 0; j < N / PP; j++)
		{
			C[i][j] = 0;

		}
    }

    pra = AA;
    prb = BB;
    psa = A;
    psb = B;

	if (id_procesu == 0) 
    {
        startwtime2 = MPI_Wtime();//czas w sekundach
    }

    for (int glowna = 0; glowna < PP; glowna++)
    {
        for (int i = 0; i < N / PP; i++)
        {   
            for (int k = 0; k < N / PP; k++)
            {
                for (int j = 0; j < N / PP; j++)
                {
                    C[i][j] += psa[i][k] * psb[k][j];
                }
            }
        }
    
        //KOMUNIKAJCA
        // MPI_Irecv(adres, ile_słów, typ_danych, odbiorca/nadawca, znacznik, zakres_procesów, Id_komunikacji);
        MPI_Irecv(pra, N * N / P, MPI_FLOAT, Prawo, tag, MPI_COMM_WORLD, &reqRecv[0]);
        MPI_Irecv(prb, N * N / P, MPI_FLOAT, Dol, tag, MPI_COMM_WORLD, &reqRecv[1]);
        MPI_Isend(psa, N * N / P, MPI_FLOAT, Lewo, tag, MPI_COMM_WORLD, &reqSend[0]);
        MPI_Isend(psb, N * N / P, MPI_FLOAT, Gora, tag, MPI_COMM_WORLD, &reqSend[1]);

        //OCZEKIWANIE NA KOMUNIKACJĘ ASYNCHRONICZNĄ
        MPI_Waitall(2, reqSend, statSend);
        MPI_Waitall(2, reqRecv, statRecv);

        //ZMIANA OBSZARÓW DANYCH LICZONYCH I PRZESYAŁANYCH
        if (mod = ((mod + 1) % 2))
        { 
            pra = A;
            prb = B;
            psa = AA;
            psb = BB;
        }
        else 
        { 
            pra = AA;
            prb = BB;
            psa = A;
            psb = B;
        }
    }
    

	if (id_procesu == 0)
	{
		endwtime = MPI_Wtime();
        printf("--------------------------------------------------------------------------------\n");
		printf("Calkowity czas przetwarzania rownoleglego wynosi: %f sekund\n", endwtime - startwtime1);
		printf("Calkowity czas obliczen rownoleglych wynosi:      %f sekund\n", endwtime - startwtime2);
        printf("--------------------------------------------------------------------------------\n");
	}

    
    MPI_Isend(C, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqSendZbierz[id_procesu]);

	if (id_procesu == 0)
	{

        // ODEBRANIE WYNIKÓW OD POZOSTAŁYCH PROCESÓW W PROCESIE 0
		for (int i = 0; i < P; i++)
        {
		   
            MPI_Irecv(C_Temp, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqRecvZbierz[i]);
            
            MPI_Wait(&reqRecvZbierz[i], &statRecvZbierz[i]);
           
            
            for (int j = 0; j < N / PP; j++)
            {
                for (int k = 0; k < N / PP; k++)
                {
                    C_Row[(i / PP) * (N / PP) + j][(i % PP) * (N / PP) + k] = C_Temp[j][k];
                }
            }
        }    

       // printf("Odebrano wyniki od wszystkich procesów\n");
        // przygotowanie C_Sek
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                C_Sek[i][j] = 0;
            }
        }


        startwtime3 = MPI_Wtime();

        // SEKWENCYJNE OBLICZENIA
		for (int i = 0; i < N; i++) 
        {
            for (int k = 0; k < N; k++) 
            {
                for (int j = 0; j < N; j++) 
                {
                    C_Sek[i][j] += A_Odczyt[i][k] * B_Odczyt[k][j];
                }
            }
        }

        endwtime = MPI_Wtime();
        printf("Calkowity czas obliczen sekwencyjnych wynosi:     %f sekund\n", endwtime - startwtime3);
        printf("--------------------------------------------------------------------------------\n");

		// porównanie poprawności obliczeń (Csek, Cglob) przy uwzględniniu progu poprawności
        int bledy = 0;
        for (int i = 0; i < N; i++)
        	{
            		for (int j = 0; j < N; j++)
            		{
                		if (C_Sek[i][j] != C_Row[i][j] && (C_Sek[i][j] / C_Row[i][j] < 0.9 || C_Sek[i][j] / C_Row[i][j] > 1.1))
                		{
                    			bledy++;
                		}
            		}
        	}
        printf("Liczba błędów:  %d\n", bledy);
        printf("Procent błędów: %.2f %%\n", (float) bledy / (N * N) * 100);
        printf("--------------------------------------------------------------------------------\n");
        // zapisanie wyników C_Row do pliku
        plik_out = fopen("Rownolegle_wynik.txt", "w");
        if (plik_out == NULL)
        {
            printf("Blad otwarcia pliku \"Rownolegle_wynik.txt\"\n");
        }
        else
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    fprintf(plik_out, "%.1f ", C_Row[i][j]);
                }
                fprintf(plik_out, "\n");
            }

            fclose(plik_out);
            printf("Wyniki zapisano do pliku \"Rownolegle_wynik.txt\"\n");
        }

        // zapisanie wyników C_Sek do pliku
        plik_out = fopen("Sekwencyjne_wynik.txt", "w"); 
        if (plik_out == NULL)
        {
            printf("Blad otwarcia pliku \"Sekwencyjne_wynik.txt\"\n");
        }
        else
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    fprintf(plik_out, "%.1f ", C_Sek[i][j]);
                }
                fprintf(plik_out, "\n");
            }

            fclose(plik_out);
            printf("Wyniki zapisano do pliku \"Sekwencyjne_wynik.txt\"\n");
            printf("--------------------------------------------------------------------------------\n");
        }
	}

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
