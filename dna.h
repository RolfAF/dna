#ifndef DNA_H
#define DNA_H

#include "utilidades.h"

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <omp.h>
#include <math.h>
#include <time.h>
#include <random>

using std::string;
using std::cout;
using std::endl;
using std::map;

char nucleotideos[] = {'A', 'C', 'G', 'T'};

string aleatorio(long tamanho, char* nucleotideos){
    cout << "Gerando sequencia aleatoria..." << endl;
    srand(time(NULL));
    string sequencia_aleatoria;
    char array[10];
    cout << array << endl;
    for(long i = 0; i < tamanho; i++){
        int numero_aleatorio = rand()%4;
        sequencia_aleatoria += nucleotideos[numero_aleatorio];
    }
    return sequencia_aleatoria;
}

string aleatorio_2(int tamanho, char* nucleotideos){
    cout << "Gerando sequencia aleatoria..." << endl;
    std::random_device dev;
    std::mt19937 aleatorio(dev());
    std::uniform_int_distribution<std::mt19937::result_type> numero(0,3);
    string sequencia_aleatoria;
    char array[10];
    cout << array << endl;
    for(int i = 0; i < tamanho; i++){
        int numero_aleatorio = numero(aleatorio);
        sequencia_aleatoria += nucleotideos[numero_aleatorio];
    }
    return sequencia_aleatoria;
}

string aleatorio_multi(int tamanho, char* nucleotideos){
    srand(time(NULL));
    cout << "Gerando sequencia aleatoria..." << endl;
    string sequencia_aleatoria[numero_threads];
    int numero_aleatorio;
    #pragma omp parallel shared(numero_aleatorio)                                               // cada thread executa uma vez
    {
        char* nuc = nucleotideos;
        int fragmento = tamanho / numero_threads;    // divide a sequencia de rna para as quatro threads
        int thread = omp_get_thread_num();                                  // numero da thread local
        int i = fragmento * thread;                                 // cada thread inicia em uma posicao diferente
        int fim_fragmento = i + fragmento;
        for(int j = i; j < fim_fragmento; j++){
            #pragma omp critical
            {
            numero_aleatorio = rand()%4;
            }
            sequencia_aleatoria[thread] += nuc[numero_aleatorio];
        }
    }
    for(int i = 0; i < 4; i++){
        cout << sequencia_aleatoria[i].length() << endl;
    }
    string sequencia = sequencia_aleatoria[0] + sequencia_aleatoria[1] + sequencia_aleatoria[2] + sequencia_aleatoria[3];
    cout << sequencia.length();
    return sequencia;
}

string valida(string sequencia_dna, char* nucleotideos){
    cout << "Validando sequencia..." << endl;
    bool sequencia_valida = true;
    #pragma omp parallel for
    for(long i = 0; i < sequencia_dna.length(); i++){
        sequencia_dna[i] = toupper(sequencia_dna[i]);
        bool caracter_valido = false;
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                caracter_valido = true;
                j = 4;
            }
        }
        if(!caracter_valido){
            cout << vermelho << "Erro, caracter invalido na sequencia de DNA: " << sequencia_dna[i] << reset << endl;
            i = sequencia_dna.length();
            sequencia_valida = false;
        }
    }
    if(sequencia_valida){
        cout << verde << "Sequencia de DNA valida ✓" << reset << endl << endl;
        return sequencia_dna;
    }else{
        return 0;
    }
}

string valida_single(string sequencia_dna, char* nucleotideos){
    cout << "Validando sequencia..." << endl;
    bool sequencia_valida = true;
    for(long i = 0; i < sequencia_dna.length(); i++){
        sequencia_dna[i] = toupper(sequencia_dna[i]);
        bool caracter_valido = false;
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                caracter_valido = true;
                j = 4;
            }
        }
        if(!caracter_valido){
            cout << vermelho << "Erro, caracter invalido na sequencia de DNA: " << sequencia_dna[i] << reset << endl;
            i = sequencia_dna.length();
            sequencia_valida = false;
        }
    }
    if(sequencia_valida){
        cout << verde << "Sequencia de DNA valida ✓" << reset << endl << endl;
        return sequencia_dna;
    }else{
        return 0;
    }
}

int * frequencia_nucleotideos(string sequencia_dna, char* nucleotideos){
    int frequencia[4][4] = {0};
    #pragma omp parallel for //private(nucleotideo)
    for(long i = 0; i < sequencia_dna.length(); i++){        // checa todas as letras da sequencia de dna
        int thread = omp_get_thread_num();
        switch(sequencia_dna[i]){
            case 'A': frequencia[thread][0]++; break;
            case 'C': frequencia[thread][1]++; break;
            case 'G': frequencia[thread][2]++; break;
            case 'T': frequencia[thread][3]++; break;
        }
    }
    static int total_frequencia[] = {0, 0, 0, 0};
    for(int i = 0; i < numero_threads; i++){                // soma os valores dos quatro arrays
        for(int j = 0; j < numero_threads; j++){
            total_frequencia[i] += frequencia[j][i];
        }
    }
    
    return total_frequencia;
}

int * frequencia_nucleotideos_2(string sequencia_dna, char* nucleotideos){
    static int frequencia[4] = {0};
    #pragma omp parallel for shared(frequencia)
    for(int i = 0; i < sequencia_dna.length(); i++){        // checa todas as letras da sequencia de dna
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                frequencia[j]++;
                j = 4;
            }    
        }
    }
    return frequencia;
}

int * frequencia_nucleotideos_single(string sequencia_dna, char* nucleotideos){
    static int frequencia[4] = {0};
    for(int i = 0; i < sequencia_dna.length(); i++){        // checa todas as letras da sequencia de dna
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                frequencia[j]++;
                j = 4;
            }    
        }
    }
    return frequencia;
}

string complemento_dna(string sequencia_dna){
    string complemento_dna = sequencia_dna;
    #pragma omp parallel for
    for(long i = 0; i < sequencia_dna.length(); i++){
        switch(sequencia_dna[i]){
            case 'A':
                complemento_dna[i] = 'T';
                break;
            case 'G':
                complemento_dna[i] = 'C';
                break;
            case 'C':
                complemento_dna[i] = 'G';
                break;
            default:
                complemento_dna[i] = 'A';
        }
    }
    return complemento_dna;
}

string complemento_dna_single(string sequencia_dna){
    string complemento_dna = sequencia_dna;
    for(long i = 0; i < sequencia_dna.length(); i++){
        switch(sequencia_dna[i]){
            case 'A':
                complemento_dna[i] = 'T';
                break;
            case 'G':
                complemento_dna[i] = 'C';
                break;
            case 'C':
                complemento_dna[i] = 'G';
                break;
            default:
                complemento_dna[i] = 'A';
        }
    }
    return complemento_dna;
}

string inverso(string sequencia){
    string sequencia_inversa = sequencia;
    long tamanho = sequencia.length();
    #pragma omp parallel for
        for(long i = 0; i < tamanho; i++){
            sequencia_inversa[tamanho-i-1] = sequencia[i];
        }    
    return sequencia_inversa;
}

string inverso_single(string sequencia){
    string sequencia_inversa = sequencia;
    long tamanho = sequencia.length();
        for(long i = 0; i < tamanho; i++){
            sequencia_inversa[tamanho-i-1] = sequencia[i];
        }    
    return sequencia_inversa;
}

string traduzir_dna_rna(string sequencia_dna){
    string sequencia_rna = sequencia_dna;
    #pragma omp parallel for
    for(long i = 0; i < sequencia_dna.length(); i++){ //checa todas as letras da sequencia de dna
        if(sequencia_dna[i] == 'T'){
            sequencia_rna[i] = 'U';
        }
    }
    return sequencia_rna;
}

string traduzir_dna_rna_single(string sequencia_dna){
    string sequencia_rna = sequencia_dna;
    for(long i = 0; i < sequencia_dna.length(); i++){ //checa todas as letras da sequencia de dna
        if(sequencia_dna[i] == 'T'){
            sequencia_rna[i] = 'U';
        }
    }
    return sequencia_rna;
}

void traduzir_rna_proteina(string sequencia_rna){
    map<string, string> codons;
        codons["UUU"] = "F";    codons["UCU"] = "S";    codons["UAU"] = "Y";    codons["UGU"] = "C";
        codons["UUC"] = "F";    codons["UCC"] = "S";    codons["UAC"] = "Y";    codons["UGC"] = "C";
        codons["UUA"] = "L";    codons["UCA"] = "S";    codons["UAA"] = "/";    codons["UGA"] = "/";
        codons["UUG"] = "L";    codons["UCG"] = "S";    codons["UAG"] = "/";    codons["UGG"] = "W";
        
        codons["CUU"] = "L";    codons["CCU"] = "P";    codons["CAU"] = "H";    codons["CGU"] = "R";
        codons["CUC"] = "L";    codons["CCC"] = "P";    codons["CAC"] = "H";    codons["CGC"] = "R";
        codons["CUA"] = "L";    codons["CCA"] = "P";    codons["CAA"] = "Q";    codons["CGA"] = "R";
        codons["CUG"] = "L";    codons["CCG"] = "P";    codons["CAG"] = "Q";    codons["CGG"] = "R";
                                                             
        codons["AUU"] = "I";    codons["ACU"] = "T";    codons["AAU"] = "N";    codons["AGU"] = "S";
        codons["AUC"] = "I";    codons["ACC"] = "T";    codons["AAC"] = "N";    codons["AGC"] = "S";
        codons["AUA"] = "I";    codons["ACA"] = "T";    codons["AAA"] = "K";    codons["AGA"] = "R";
        codons["AUG"] = "M";    codons["ACG"] = "T";    codons["AAG"] = "K";    codons["AGG"] = "R";
                                                                                   
        codons["GUU"] = "V";    codons["GCU"] = "A";    codons["GAU"] = "D";    codons["GGU"] = "G";
        codons["GUC"] = "V";    codons["GCC"] = "A";    codons["GAC"] = "D";    codons["GGC"] = "G";
        codons["GUA"] = "V";    codons["GCA"] = "A";    codons["GAA"] = "E";    codons["GGA"] = "G";
        codons["GUG"] = "V";    codons["GCG"] = "A";    codons["GAG"] = "E";    codons["GGG"] = "G";
    string proteinas[numero_threads];
    #pragma omp parallel                                                // cada thread executa uma vez
    {
        long thread = omp_get_thread_num();                                          // numero da thread local
        long tamanho_fragmento = sequencia_rna.length() / numero_threads;            // divide a sequencia de rna para as quatro threads
        long i = tamanho_fragmento * thread;                                         // cada thread inicia em uma posicao diferente
        while(i%3 != 0){                                                            // mantem sincronizacao (reading frame)
            i++;
        }
        long fim_fragmento;
        if(thread == 3){
            tamanho_fragmento = sequencia_rna.length() - (tamanho_fragmento * 3);   // ultima thread recebe o resto da divisao
            fim_fragmento = sequencia_rna.length();
        }else{
            fim_fragmento = i + tamanho_fragmento - 1;
        }
        bool traduzindo_proteina;
        string proteina;
        while(i+2 < fim_fragmento || (traduzindo_proteina && thread != 3)){                  // while pode continuar caso proteina continue alem do fim do fragmento
            string codon = sequencia_rna.substr(i,3);
            //cout << colorir(codon) << " ";
            if(codon == "AUG"){                                       // codon de inicio
                traduzindo_proteina = true;
            }else if(traduzindo_proteina && (codon == "UAA" || codon == "UAG" || codon == "UGA")){          // codon de fim
                traduzindo_proteina = false;
                proteina += " ";                                            // espaco vazio delimita fim da proteina
                proteinas[thread] += proteina;                      // vetor de proteinas recebe nova proteina
                proteina = "";                                              // limpa a variavel para receber nova proteina
            }
            if(traduzindo_proteina){
                proteina += codons[codon];                                  // amino-acido correspondente adicionado a proteina
            }
            i+=3;
        }
    }
    cout << endl;
    for(long i = 0; i < numero_threads; i++){
        cout << i << ": " << endl;
        cout << proteinas[i] << " " << endl;
    }
    
}

void traduzir_rna_proteina_branchless(string sequencia_rna){
    string proteinas[numero_threads];
    #pragma omp parallel                                                // cada thread executa uma vez
    {
        long thread = omp_get_thread_num();                                          // numero da thread local
        long tamanho_fragmento = sequencia_rna.length() / numero_threads;            // divide a sequencia de rna para as quatro threads
        long i = tamanho_fragmento * thread;                                         // cada thread inicia em uma posicao diferente
        while(i%3 != 0){                                                            // mantem sincronizacao (reading frame)
            i++;
        }
        long fim_fragmento;
        if(thread == 3){
            tamanho_fragmento = sequencia_rna.length() - (tamanho_fragmento * 3);   // ultima thread recebe o resto da divisao
            fim_fragmento = sequencia_rna.length();
        }else{
            fim_fragmento = i + tamanho_fragmento - 1;
        }
        bool traduzindo_proteina;
        string proteina;
        while(i+2 < fim_fragmento || (traduzindo_proteina && thread != 3)){                  // while pode continuar caso proteina continue alem do fim do fragmento
            string codon = sequencia_rna.substr(i,3);
            //cout << colorir(codon) << " ";
            if(codon == "AUG"){                                       // codon de inicio
                traduzindo_proteina = true;
            }
            if(traduzindo_proteina){   
                //F 70 L 76 I 73 M 77 V 86 S 83 P 80 T 84 A 65 Y 89 H 72 Q 81 N 78 K 75 D 68 E 69 C 67 W 87 R 82 S 83 G 71 / 47 ' ' 32
                proteina += (char) ( (codon == "UUU" || codon == "UUC")*70    
                                   + (codon == "UUA" || codon == "UUG")*76    
                                   + (codon == "UCU" || codon == "UCC" || codon == "UCA" || codon == "UCG")*83
                                   + (codon == "UAU" || codon == "UAC")*89    
                                   + (codon == "UGU" || codon == "UGC")*67                              
                                   + (codon == "UGG")*87                                 
                                   + (codon == "CUU" || codon == "CUC" || codon == "CUA" || codon == "CUG")*76    
                                   + (codon == "CCU" || codon == "CCC" || codon == "CCA" || codon == "CCG")*80    
                                   + (codon == "CAU" || codon == "CAC")*72    
                                   + (codon == "CAA" || codon == "CAG")*81                                                                                                         
                                   + (codon == "CGU" || codon == "CGC" || codon == "CGA" || codon == "CGG")*82                                                                                                         
                                   + (codon == "AUU" || codon == "AUC" || codon == "AUA")*73  
                                   + (codon == "AUG")*77
                                   + (codon == "ACU" || codon == "ACC" || codon == "ACA" || codon == "ACG")*84    
                                   + (codon == "AAU" || codon == "AAC")*78                                                                               
                                   + (codon == "AAA" || codon == "AAG")*75      
                                   + (codon == "AGU" || codon == "AGC")*83
                                   + (codon == "AGA" || codon == "AGG")*82
                                   + (codon == "GUU" || codon == "GUC" || codon == "GUA" || codon == "GUG")*86    
                                   + (codon == "GCU" || codon == "GCC" || codon == "GCA" || codon == "GCG")*65    
                                   + (codon == "GAU" || codon == "GAC")*68    
                                   + (codon == "GAA" || codon == "GAG")*69        
                                   + (codon == "GGU" || codon == "GGC" || codon == "GGA" || codon == "GGG")*71 );             
                if(codon == "UAA" || codon == "UAG" || codon == "UGA"){
                    traduzindo_proteina = false;
                    proteina += " ";
                    proteinas[thread] += proteina;
                    proteina = "";
                }                
            }
            i+=3;
        }
    }
    cout << endl;
    for(long i = 0; i < numero_threads; i++){
        cout << i << ": " << endl;
        cout << proteinas[i] << " " << endl;
    }
    
}

void traduzir_rna_proteina_single(string sequencia_rna){
    std::vector<string> proteinas;
    string proteina;
    long encontrou_proteina = 0;
    map<string, string> codons;
        codons["UUU"] = "F";    codons["UCU"] = "S";    codons["UAU"] = "Y";    codons["UGU"] = "C";
        codons["UUC"] = "F";    codons["UCC"] = "S";    codons["UAC"] = "Y";    codons["UGC"] = "C";
        codons["UUA"] = "L";    codons["UCA"] = "S";    codons["UAA"] = "/";    codons["UGA"] = "/";
        codons["UUG"] = "L";    codons["UCG"] = "S";    codons["UAG"] = "/";    codons["UGG"] = "W";
        
        codons["CUU"] = "L";    codons["CCU"] = "P";    codons["CAU"] = "H";    codons["CGU"] = "R";
        codons["CUC"] = "L";    codons["CCC"] = "P";    codons["CAC"] = "H";    codons["CGC"] = "R";
        codons["CUA"] = "L";    codons["CCA"] = "P";    codons["CAA"] = "Q";    codons["CGA"] = "R";
        codons["CUG"] = "L";    codons["CCG"] = "P";    codons["CAG"] = "Q";    codons["CGG"] = "R";
                                                             
        codons["AUU"] = "I";    codons["ACU"] = "T";    codons["AAU"] = "N";    codons["AGU"] = "S";
        codons["AUC"] = "I";    codons["ACC"] = "T";    codons["AAC"] = "N";    codons["AGC"] = "S";
        codons["AUA"] = "I";    codons["ACA"] = "T";    codons["AAA"] = "K";    codons["AGA"] = "R";
        codons["AUG"] = "M";    codons["ACG"] = "T";    codons["AAG"] = "K";    codons["AGG"] = "R";
                                                                                   
        codons["GUU"] = "V";    codons["GCU"] = "A";    codons["GAU"] = "D";    codons["GGU"] = "G";
        codons["GUC"] = "V";    codons["GCC"] = "A";    codons["GAC"] = "D";    codons["GGC"] = "G";
        codons["GUA"] = "V";    codons["GCA"] = "A";    codons["GAA"] = "E";    codons["GGA"] = "G";
        codons["GUG"] = "V";    codons["GCG"] = "A";    codons["GAG"] = "E";    codons["GGG"] = "G";
    for(long i = 0; i+2 < sequencia_rna.length(); i+=3){
        string codon = sequencia_rna.substr(i,3);
        if(codon == "AUG"){                                   // codon de inicio
            encontrou_proteina = 1;
        }else if(encontrou_proteina == 1 && (codon == "UAA" || codon == "UAG" || codon == "UGA")){  // codon de fim
            encontrou_proteina = 0;
            proteinas.push_back(proteina);                          // vetor de proteinas recebe nova proteina
            proteina = "";                                          // limpa a variavel para receber nova proteina
        }
        if(encontrou_proteina == 1){
            proteina += codons[codon];                              // amino-acido correspondente adicionado a proteina
        }        
    }
    cout << endl;
    for(long i = 0; i < proteinas.size(); i++){
        cout << proteinas[i] << " ";
    }
    cout << endl << endl;
}

void gerar_codons(string sequencia_dna){
    string codons[6];
    string sequencia_dna_inversa = inverso(sequencia_dna);
    codons[0] = sequencia_dna;
    codons[1] = sequencia_dna.substr(1,sequencia_dna.length());
    codons[2] = sequencia_dna.substr(2,sequencia_dna.length());
    codons[3] = sequencia_dna_inversa;
    codons[4] = sequencia_dna_inversa.substr(2,sequencia_dna_inversa.length());
    codons[5] = sequencia_dna_inversa.substr(3,sequencia_dna_inversa.length());
    cout << "codon 1: " << codons[0] << endl;
    cout << "codon 2: " << codons[1] << endl;
    cout << "codon 3: " << codons[2] << endl;
    cout << "codon 4: " << codons[3] << endl;
    cout << "codon 5: " << codons[4] << endl;
    cout << "codon 6: " << codons[5] << endl;
    
}

void encontra_assinatura(string sequencia, string assinatura){
    long posicao[sequencia.length()];
    cout << sequencia << endl;
    //#pragma omp parallel for
    for(long i = 0; i < sequencia.length(); i++){
        for(long j = 0; j < assinatura.length(); j++){
            if(assinatura[j] != sequencia[i+j]){
                break;
            }
            if(j == assinatura.length()-1){
                //#pragma omp critical
                //cout << i+1 << " ";
                posicao[i] = 1;
            }
        }
    }
    for(long i = 0; i < sequencia.length(); i++){   
        if(posicao[i] == 1){
           cout << i+1 << " ";
        }
    }
    cout << endl << endl;
}

#endif
