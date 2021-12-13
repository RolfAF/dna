#ifndef UTILIDADES_H
#define UTILIDADES_H

//#include "dna.h"

#include <iostream>
#include <string>

using std::string;
using std::cout;
using std::endl;

string vermelho="\e[38;2;255;100;100m";     string verde="\e[38;2;53;255;145m";     string azul="\e[38;2;175;175;255m";
string amarelo="\e[38;2;255;225;0m";        string reset="\e[0m";                   string bold="\e[1m";
int tamanho_linha = 100;
int numero_threads = 4;

string colorir(string sequencia){
    int tamanho;
    string sequencia_colorida;
    //sequencia_colorida = bold;
    string reticencia = "";
    if(sequencia.length() <= tamanho_linha){
        tamanho = sequencia.length();
    }else{
        tamanho = tamanho_linha;
        reticencia = " ...";
    }    
    for(int i = 0; i < tamanho; i++){
        switch(sequencia[i]){
            case 'A':
                sequencia_colorida += verde + sequencia[i];
                break;
            case 'C':
                sequencia_colorida += azul + sequencia[i];
                break;
            case 'G':
                sequencia_colorida += amarelo + sequencia[i];
                break;
            case 'T':
                sequencia_colorida += vermelho + sequencia[i];
                break;
            case 'U':
                sequencia_colorida += vermelho + sequencia[i];
                break;
        }
    }
    sequencia_colorida += reset += reticencia;
    return sequencia_colorida;
}

void imprime_complemento_dna(string sequencia_dna, string complemento_dna){
    int tamanho = tamanho_linha;
    if(sequencia_dna.length() <= tamanho_linha){
        tamanho = sequencia_dna.length();
    }
    cout << "DNA: " << colorir(sequencia_dna.substr(0,tamanho_linha)) << endl;
    string aux;
    for(int i = 0; i < tamanho; i++){
        aux += "￬";
    }
    cout << "     " << aux << endl;
    cout << "DNA: " << colorir(complemento_dna.substr(0,tamanho_linha)) << endl << endl;
}

void imprime_traducao_dna_rna(string sequencia_dna, string sequencia_rna){
    int tamanho = tamanho_linha;
    if(sequencia_dna.length() <= tamanho_linha){
        tamanho = sequencia_dna.length();
    }
    cout << "DNA: " << colorir(sequencia_dna.substr(0,tamanho_linha)) << endl;
    string aux;
    for(int i = 0; i < tamanho; i++){
        if(sequencia_dna[i] == 'T'){
            aux += "￬";
        }else{
            aux += "|";
        }
    }
    cout << "     " << aux << endl;
    cout << "RNA: " << colorir(sequencia_rna.substr(0,tamanho_linha)) << endl << endl;
}

#endif
