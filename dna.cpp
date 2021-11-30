#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <omp.h>
#include <math.h>
#include <time.h>
//#include <algorithm> 
//#include <cstdio>
//#include <random>
//using namespace std;
using std::string;
using std::cout;
using std::endl;
using std::map;

string vermelho="\e[38;2;255;100;100m";     string verde="\e[38;2;53;255;145m";     string azul="\e[38;2;175;175;255m";
string amarelo="\e[38;2;255;225;0m";        string reset="\e[0m";                   string bold="\e[1m";
int tamanho_linha = 75;
int numero_threads = 4;

char nucleotideos[] = {'A', 'C', 'G', 'T'};

string aleatorio(int tamanho, char* nucleotideos){
    cout << "Gerando sequencia aleatoria..." << endl;
    int min = 0;
    int max = 3;
    srand(time(NULL));
    string sequencia_aleatoria;
    char array[10];
    cout << array << endl;
    for(int i = 0; i < tamanho; i++){
        int numero_aleatorio = rand()%(max-min+1)+min;
        sequencia_aleatoria += nucleotideos[numero_aleatorio];
    }
    return sequencia_aleatoria;
}

string aleatorio_multi(int tamanho, char* nucleotideos){
    srand(time(NULL));
    cout << "Gerando sequencia aleatoria..." << endl;
    string sequencia_aleatoria[numero_threads];
    #pragma omp parallel                                                // cada thread executa uma vez
    {
        char* nuc = nucleotideos;
        int fragmento = tamanho / numero_threads;    // divide a sequencia de rna para as quatro threads
        int thread = omp_get_thread_num();                                  // numero da thread local
        int i = fragmento * thread;                                 // cada thread inicia em uma posicao diferente
        int fim_fragmento = i + fragmento;
        int min = 0;
        int max = 3;
        for(int j = i; j < fim_fragmento; j++){
            int numero_aleatorio = rand()%(max-min+1)+min;
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
    #pragma omp parallel for
    for(int i = 0; i < sequencia_dna.length(); i++){
        sequencia_dna[i] = toupper(sequencia_dna[i]);
        bool caracter_valido = false;
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                caracter_valido = true;
            }
        }
        if(!caracter_valido){
            cout << vermelho << "Erro, caracter invalido na sequencia de DNA: " << sequencia_dna[i] << reset << endl;
            i = sequencia_dna.length();
        }
    }
    cout << verde << "Sequencia de DNA valida ✓" << reset << endl << endl;
    return sequencia_dna;
}

int * frequencia_nucleotideos(string sequencia_dna, char* nucleotideos){
    int frequencia[4][4] = {0};
    #pragma omp parallel for
    for(int i = 0; i < sequencia_dna.length(); i++){        // checa todas as letras da sequencia de dna
        int thread = omp_get_thread_num();
        for(int j = 0; j < 4; j++){                         // compara a letra com os quatro possiveis nucleotideos
            if(sequencia_dna[i] == nucleotideos[j]){
                frequencia[thread][j]++;
                //break;
                j = 4;
            }    
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

string traduzir_dna_rna(string sequencia_dna){
    string sequencia_rna = sequencia_dna;
    map<char, char> complemento;
        complemento['A'] = 'A';    complemento['C'] = 'C';
        complemento['G'] = 'G';    complemento['T'] = 'U';
    #pragma omp parallel for
    for(int i = 0; i < sequencia_dna.length(); i++){ //checa todas as letras da sequencia de dna
        sequencia_rna[i] = complemento[sequencia_dna[i]];
    }
    return sequencia_rna;
}

string complemento_dna(string sequencia_dna){
    string complemento_dna = sequencia_dna;
    map<char, char> complemento;
        complemento['A'] = 'T';    complemento['C'] = 'G';
        complemento['G'] = 'C';    complemento['T'] = 'A';
    #pragma omp parallel for
    for(int i = 0; i < sequencia_dna.length(); i++){
        complemento_dna[i] = complemento[sequencia_dna[i]];
    }
    return complemento_dna;
}

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

string inverso(string sequencia){
    string sequencia_inversa[numero_threads];
    #pragma omp parallel                                                // cada thread executa uma vez
    {
        int thread = omp_get_thread_num();                                          // numero da thread local
        int tamanho_fragmento = sequencia.length() / numero_threads;            // divide a sequencia para as quatro threads
        int inicio_fragmento = tamanho_fragmento * thread;                                         // cada thread inicia em uma posicao diferente
        int fim_fragmento;
        if(thread == 3){
            tamanho_fragmento = sequencia.length() - (tamanho_fragmento * 3);   // ultima thread recebe o resto da divisao
            fim_fragmento = sequencia.length();
        }else{
            fim_fragmento = inicio_fragmento + tamanho_fragmento - 1;
        }
        for(int j = fim_fragmento; j >= inicio_fragmento; j--){
            sequencia_inversa[thread] += sequencia[j];
        }
    }
    return sequencia_inversa[3] + sequencia_inversa[2] + sequencia_inversa[1] + sequencia_inversa[0];
}

void imprime_complemento_dna(string sequencia_dna){
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
    cout << "DNA: " << colorir(complemento_dna(sequencia_dna.substr(0,tamanho_linha))) << endl << endl;
}

void imprime_traducao_dna_rna(string sequencia_dna){
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
    cout << "RNA: " << colorir(traduzir_dna_rna(sequencia_dna.substr(0,tamanho_linha))) << endl << endl;
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
        int thread = omp_get_thread_num();                                          // numero da thread local
        int tamanho_fragmento = sequencia_rna.length() / numero_threads;            // divide a sequencia de rna para as quatro threads
        int i = tamanho_fragmento * thread;                                         // cada thread inicia em uma posicao diferente
        while(i%3 != 0){                                                            // mantem sincronizacao (reading frame)
            i++;
        }
        int fim_fragmento;
        if(thread == 3){
            tamanho_fragmento = sequencia_rna.length() - (tamanho_fragmento * 3);   // ultima thread recebe o resto da divisao
            fim_fragmento = sequencia_rna.length();
        }else{
            fim_fragmento = i + tamanho_fragmento - 1;
        }
        bool traduzindo_proteina;
        string proteina;
        while(i+2 < fim_fragmento || traduzindo_proteina){                  // while pode continuar caso proteina continue alem do fim do fragmento
            string codon = sequencia_rna.substr(i,3);
            //cout << colorir(codon) << " ";
            if(codons[codon] == "M"){                                       // codon de inicio
                traduzindo_proteina = true;
            }else if(traduzindo_proteina && codons[codon] == "/"){          // codon de fim
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
    for(int i = 0; i < numero_threads; i++){
        cout << i << ": " << endl;
        cout << proteinas[i] << " " << endl;
    }
    
}

void traduzir_rna_proteina_single(string sequencia_rna){
    std::vector<string> proteinas;
    string proteina;
    int encontrou_proteina = 0;
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
    for(int i = 0; i+2 < sequencia_rna.length(); i+=3){
        string codon = sequencia_rna.substr(i,3);
        if(codons[codon] == "M"){                                   // codon de inicio
            encontrou_proteina = 1;
        }else if(encontrou_proteina == 1 && codons[codon] == "/"){  // codon de fim
            encontrou_proteina = 0;
            proteinas.push_back(proteina);                          // vetor de proteinas recebe nova proteina
            proteina = "";                                          // limpa a variavel para receber nova proteina
        }
        if(encontrou_proteina == 1){
            proteina += codons[codon];                              // amino-acido correspondente adicionado a proteina
        }        
    }
    cout << endl;
    for(int i = 0; i < proteinas.size(); i++){
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


int main(int argc, char *argv[]){
    
    
    //string sequencia_dna = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"; // rosalind 1
    //string sequencia_dna = "ATCCACCCTCCGTCAGCTCTGCCACTCAGGGCAATAGAAAAATGTAGGAAGGGTACCCGAAGGCCGTTGCATTGGATGTTTGTGATAGCGTGATGGATGTCGGGTGGATAGATAAGCCCCGCCACTCGGTGGATATCAGCGAACACAGTAGCAACTGTTTGTGACATCTACCAACAATAGTCATTTTCGCCATTAATTTACGCACCCAGGTGGATCCCAACTGGTCACAAGTTGTTCTGCGTACGAATCGTAACGGGAACGGGTTCGTATTTATTTTTTGGACCGATCCCACGTAGATACGGTGGCCTCCTAATGCTACTATATATGCGCCGTATCGTGGAATCCAAACCGCCCAGATAAAGCGACGGATCCGCCTATTTATAGCGGAGTCGATTAAAGGTAAAGAGCAATCAATACAAGACTGAACCATCAGTTCCCAGGAGAAGAGCAACGTGTACACGGCTTCAGAGTGAGCTACACTGGATTGTACGTCCAGCTGTATCAGCCATTCTGTCGCGGGTGTCCGTCTATATTTAACGAATCGCAATCATTTGAGAATCTGGCCCTACGAAGATGCGATTCGATAATCAGATGTGTTTGTCATGGGACGGAATCCTATCAACGGTTTGTGATGCTCAGAAACATGCTATTGAGCATCGAAGCGTAATCCTGAAAGGGTTTAAGGGTACCATGCTCCCGACCCAAGACGTAGAAGAAGTTCAAGCAGCCTATGACTTTGTCGCGAAGATGGATACCGTAGTTGCACGCAATTTCGCGAAACTGTTGTCTCTGTTATACAAGTTGGAGGTAAATATGGCAAAACTGGTGGCTTCATGAGGATCAGAACCTCCGAGAGACGCTTAAGCTGTAATTATCCCGTCCGAGAGTCCCCTATATAAGCCCCTCCCTATCATGGATGAGCTCTGTCACAGACACCATTAGCAAACCTTGGTATGCGCTGCCTATGTTTGTATTAGGAGAGATACACCAGAGCAGCTTT"; // rosalind 2
    
    //string sequencia_dna = "ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"; // rosalind 3 proteina
    
    //string sequencia_dna = "AGAGCAGCGTGGTCCCAGGAGGCATACACAGAGAGCCACGGCCAGGGCTGAAACAGTCTGTTGAGTGCAGCCATGGGGGACGTCCTGGAACAGTTCTTCATCCTCACAGGGCTGCTGGTGTGCCTGGCCTGCCTGGCGAAGTGCGTGAGATTCTCCAGATGTGTTTTACTGAACTACTGGAAAGTTTTGCCAAAGTCTTTCTTGCGGTCAATGGGACAGTGGGCAGTGATCACTGGAGCAGGCGATGGAATTGGGAAAGCGTACTCGTTCGAGCTAGCAAAACGTGGACTCAATGTTGTCCTTATTAGCCGGACGCTGGAAAAACTAGAGGCCATTGCCACAGAGATCGAGCGGACTACAGGGAGGAGTGTGAAGATTATACAAGCAGATTTTACAAAAGATGACATCTACGAGCATATTAAAGAAAAACTTGCAGGCTTAGAAATTGGAATTTTAGTCAACAATGTCGGAATGCTTCCAAACCTTCTCCCAAGCCATTTCCTGAACGCACCGGATGAAATCCAGAGCCTCATCCATTGTAACATCACCTCCGTAGTCAAGATGACACAGCTAATTCTGAAACATATGGAATCAAGGCAGAAAGGTCTCATCCTGAACATTTCTTCTGGGATAGCCCTGTTTCCTTGGCCTCTCTACTCCATGTACTCAGCTTCCAAGGCGTTTGTGTGCGCATTTTCCAAGGCCCTGCAAGAGGAATATAAAGCAAAAGAAGTCATCATCCAGGTGCTGACCCCATATGCTGTCTCGACTGCAATGACAAAGTATCTAAATACAAATGTGATAACCAAGACTGCTGATGAGTTTGTCAAAGAGTCATTGAATTATGTCACAATTGGAGGTGAAACCTGTGGCTGCCTTGCCCATGAAATCTTGGCGGGCTTTCTGAGCCTGATCCCGGCCTGGGCCTTCTACAGCGGTGCCTTCCAAAGGCTGCTCCTGACACACTATGTGGCATACCTGAAGCTCAACACCAAGGTCAGGTAGCCAGGCGGTGAGGAGTCCAGCACAACCTTTTCCTCACCAGTCCCATGCTGGCTGAAGAGGACCAGAGGAGCAGACCAGCACTTCAACCTAGTCCGCTGAAGATGGAGGGGGCTGGGGTCACAGAGGCATAGAATACACATTTTTTGCCACTTTA"; // testosterona
    
    //MGDVLEQFFILTGLLVCLACLAKCVRFSRCVLLNYWKVLPKSFLRSMGQWAVITGAGDGIGKAYSFELAKRGLNVVLISRTLEKLEAIATEIERTTGRSVKIIQADFTKDDIYEHIKEKLAGLEIGILVNNVGMLPNLLPSHFLNAPDEIQSLIHCNITSVVKMTQLILKHMESRQKGLILNISSGIALFPWPLYSMYSASKAFVCAFSKALQEEYKAKEVIIQVLTPYAVSTAMTKYLNTNVITKTADEFVKESLNYVTIGGETCGCLAHEILAGFLSLIPAWAFYSGAFQRLLLTHYVAYLKLNTKVR // <- proteina da testosterona
    
    string sequencia_dna = aleatorio(1000000000, nucleotideos);
    
    
    sequencia_dna = valida(sequencia_dna, nucleotideos);
    cout << "Sequencia de DNA: " << colorir(sequencia_dna) << endl << endl;
    
    cout << "Frequencia de nucleotideos:" << endl;
    int *frequencia = frequencia_nucleotideos(sequencia_dna, nucleotideos);
    for(int i = 0; i < 4; i++){
        cout << nucleotideos[i] << ": " << frequencia[i] << endl;
    }
    cout << "Rosalind: ";
    for(int i = 0; i < 4; i++){
        cout << frequencia[i] << " ";
    }
    cout << endl << endl;
    
    cout << "Complemento DNA: " << endl;
    string complemento = complemento_dna(sequencia_dna);
    imprime_complemento_dna(sequencia_dna);
    cout << "Complemento reverso DNA:" << endl;
    cout << colorir(inverso(complemento)) << endl << endl;
    
    cout << "Traducao de DNA para RNA:" << endl;
    string sequencia_rna = traduzir_dna_rna(sequencia_dna);
    imprime_traducao_dna_rna(sequencia_dna);
    
    cout << "Sequencia de RNA traduzida para Proteinas:" << endl;
    traduzir_rna_proteina(sequencia_rna);
    
    //cout << "Encontra assinatura:" << endl;
    //encontra_assinatura("GATATATGCATATACTT", "ATAT");
    
    string proteina="AUGCGGAUACCGGGUGCUGUUACGUUCGUUAGGAGGUACAUAGGCAAUAUGAACAGAAGUGGUAGGCACCUUCUGUCCGUACUAUACGUUUGCUCCGGCGGGAAAAGGAUCCCCCCGCGUUACAUGCACCUGAUGGAGCUGGGAGGGACACCCGAAGUUCGUACAGAAGUAAACAGAGGACGGGGUCCAUGUAGACACACUAGCAACCCCACCGUGAAUAUUUCAAGGGGACAGCGUUUAUCAGAUAGUAAGCACUGUAGGGGUUUCUACAUGCCUUCUCUGUUAUCAGUCAACACUCAUGGUUCUGCUUCGGAUAACCGACCGGCCGGAAUUUCCUUAUGGCAACUCGAGCUUGGGACCUGGUUUACGCAGUAUCUAUUCGAAAUAAAGCUGCAUAAAAACGUGAUAGCACUCCUUUCUACAUCGCGCUUAUGGAAAACACGGAACGGAAGGUUGAGUGGCAACUACCCUUCACAUCCCCGGCAUAAUUGUUCCCCUACUAGUGGGCCAGCAAGUAACCUUGUAGCUCUCCGGCUCUCGAAGGCAGGAUUGAACCCCUCUAUAGUCGUGUUCUGCCACCGCAAGGAAUUAAAACUCGCACUUAUUUAUGUCGUCACAAAUUUGAGCCGUACAGGUAAUUGGCAGAUUAUAUUAGCUCGCCCGUGGCACGGGCCCGCCCUCUACCGUACCAAAGGUCCUACUACAGUCCCAUCAUUCCAUUUCUGCAGUCCUGUAAUGGACCGCGGGGUACUGGGCCGAACGAAAGAUUCACAAGAUCGGCAGAUAGCUGCCUCCUGCCUACCUGCACCCACCGUUAUGGUUACCCGGUCAGGCCACACAGCGCUAACGUCUGGUUGGAACAUCGACGUCAUGUAUAAGAAAUUAGCAGGCCUAAAGCAUCUCGUAGCCAGUGCAUUCCCCCAUUCCCCAAGUUCGCUAUUCAGUGAUAAAAUACAGCGCCGCGUUGUCACAUCGGGAUUUUAUGCCCCGGUCUGGGGUUUUGCAUUGGUACAGGAUACAGAGCAUAGAUACUGUUCCACUGCCCGCACCACUUCAGGUCUCUGCGACAUCUGCCAACAUUCUUACGAAUCGGAGACUGAUUCCCGUCCGUGCAGCCGAACACGACAAUCGGUACAAGGAUUGGGCAGCUAUGGCCAGCGAUGUGGCCGUCCGACGGCUCGCCAAACGUUGCCUCCAGAGGGAACAGAGCCCGGACCUGUGGGCGAUCUUGAAAGUGGCUGUAUGGAAUUUAAUGUCCGUCUGCUCGUAUUGGGGCGUCACACUCUAUCUGUGUGGUGUGGAUUGACCCGGGUCGCAAACUGCAGCGAUGGGACCAGCGUAGAUCGUUUACAAGCCGGGGAGUUACUAAUCAGUGACGUGGGCCAUUAUGUGAAUUUUUGUCUAACUUAUCUGCCCGAAACCGUCCAUCCCUGCGCGCGUGUGGGCUGUUUGCCCACGUUUAACAGGAUCAUAAGGACGUUUGUACGUUUAGUAACACCCCAUUUGAUGCUGUUUGUGUUACUCACAGCGGCUCUUAUGACUCAGGCUGUCUGUUGGGAUAAGAGAUAUAGUACGAAACGCACGCUGGAGAGCCUAUCGGGGGUGCACUACUCACAUAAUUUACAGACGACAGUGCGAGCCUGCGGUAUCCGUUGGAUGUUAUGCGGUUUUCGAACCAAGUUCCAACGUAGUGCCAAUUCGCGUACGUGUGUGUCACUUAUGGUCAUUGGAAGAAUCCAAGGUAAGAAAUUACAGGCGACCGCUAGGGACAAAAUGUUCUCCGUUGAGGCACACUUUCAAGUCCUUAAGGUAACUAACCUAUCCGCGGCAAAGCUUGAACUAAUUAGUCUAGGAAAGCAAGAGGCGGAACCCCAGCAAACGCUGACACUACAUCCACACAUUGGCGCACAAUAUCGCGUGUUCUCAGAUCCGGGUGUUCGUGGAGUGCCGAGGUUCCGUCAGCAUAACCACUAUUGCGUAAGCAACUGGGAACCAGGUACCUUUAAGGAAACGUGGGAAUUCGCACAGGAACAAAAUAAUUGGUCCUUUAGCAAAGACGAGGUGGCUUGCAAGAUACCCGAGUAUCGGUGGUGUCGGAUUGAUUCACUGCGCUCACUGCUUUAUGCCCACGUCUACGCCCACAAACUCGUGCAAUCGACUAGGUUUGAACGAUUGGCGGGGGUCGGCUGCCUACAUUUAAUCGCGCUGGGACCCGUGUGGUUGACGAAAGCUCAGACAUUGCUCACUUUGGGGAAUCGCGUGUCCGACCUUAAAGGUUGCUACAUAGGGACAAAUUUGGUGUUACUGCAUCACGUUAUUUGUUCCGAGUAUGGAUGUAAGAGCAUUACGGUACAAGUGAGGUGCUUCAUCAAUAUAGCAACUAGGAUCUGCGGGGUCUAUGUGUGGAGCUUCUUGACGACUCCGUAUGGAACACCUAGUCGCCAAGACAUCUCCAGGGUAUAUGUGGUCCCUAACACAAUGCCAUUCAUAGGUAGCGGAUUUGUGCGUGGCAAAAAGGCACUACCCCCCAACUGUAAGGGGAGAUCUGGUAAGGCAACACUCCACACUUGUGAACCUUGCUGGCUCCACGAUCCGAUCCUUGUCUAUAUAGUCGGGGCGAUGAACUGGAGACGUGUUCAGAAGAGCGUCGCGCAGUUCAGUCUGAGGGUGUCGGGAAUUCCCGGGAAUUCUUUGCGAUUUAGUAACCCUACAUGCGGUGGGCCCCGGAUUUAUACACCACAAAGCUUAGACCCAGCAUGUUGCUGCGGGCGGGAACACGAAAUGGACGUAUGCUUAGUACAAAUAGAGAUUUACGGUGAAUCCAAGUCUAUUCUACCGAAGGAUCCGGCGGGCAAUGAACUCGAACACUGUGGCUUCCCGAAUUUCUCGUGCAUUAGGCUGUUUCUCUGGAACCAUAUGUGGCGGGUGGCUUGUCAAUCGCCAGUCAUUGCUUAUUUUAAGACCGUAACUACAGGGGCUGCGCUUCGUCGGAAUACCUCUCUGAGUUGCACGAGACUGCGCCCAGACGCUGGCUGUCACCUGGCGCGGGAGUACAUAGUCUGCAUUGGUGUAAAGUUGCCAGACGAGCCGUUAGUAGAGAAAUCCCUAACCUAUACAGUACCUGAUGUCGCCCGAGUAAAACUUAUAGCCACUAAGGUCCUUGUUCGUGAAUGUUUCACCUUACAUUCUCCCAAAGGAUAUGUCCUGGACGAUCGCAUAGACUUUGUUACCGACGCAACGAAUCAAGAAACUGUGUCUUAUAUGUCAAUUUUAACUUGCUCAGAAGAUGUACUAAUAAAACGUAUGUCCCGGGUACCAUCAGUCUGCUGGCUCGCAACAUUUAAAAUUUUAGCGCCGCGAACUUGGCGCGCCAUAAUACGAAGUCGGAAUGCCGUCAGCCUAUUGUUCCGAGGUGGUAUUUUAGAACCGCGCUUGUACCUCUGCACAUGUGCUUAUGAGGGCGCAGUCAUGGAAUCAGUACAUAUGAUUGGUCAUAAGCUGGUGGGGGGGAAGGGGGCUUUCUGUAAGACGAUUGCGGCCAAGGCUCGGUACAAAGAGUGUGUUCCACCCAUACAACAAGCCAAUUCUGCUAGGUACCAAAACUAUAGAGAACUCCCAAAGGUUGCUAAAUGCUCAAGAGUAGCCCUACAAUUUGGCCUUUGGUGUGAGUUCCGUUCAUUGACUGCACCUCCGAUGCUAUAUGCCAGAGCUAACCUCCUGUCGAUCGACAAUUUUCUUCGAGCCCCGUGUUUCAGUUCGACAAUCCAUGACUCCGUCCGAUUAUUUCAUUUACAAAUGGUUAUAAAGAAAAGCAGGCGAAUCGCCGUGUGUAUUCGACCUCCCCUAAGUGAACGGCGUGGAUCCCCUUGCAUUAAUUGCGAAGGUGUCCGAACAGCUUCCAGGAAAAUCAUUACACUCCACGCCAAAAUUCUGACAGCUACUAUAGCCGUGCCAGGCUACUGCCGUGCCGAUUCAGGAUCACCUGGAGAAAAUGAAUACAGCGAAACUGUGCAAGCUCACUUGGCACAGGUUGGCGCCACGUUUAGUAAGGGUCUCCCUGGUACACGUACCGCGUGCUGGGUUCCAUCGCCUAAGGUCAAUUUAACAUCGAUACCUAAGGUUCAGUAUAUAAGUUUUAGGCGGAAGUCUUGGGAGCUUCAUUGCGCUAGUCAACCGCUGGCAGUAUGGUCCGUUGAUUCCGUUCACUGGGGGCUGCAACGGUACAUGCUGCCAGUCAACAAAAAGAGUCGUAAAGUCCUCGAGGUUAACAAUAGUUCUAUAAUUACGGCACCUAACGACCAGCACUCUUACCAUGCUAUAGGUCACUGGACCACACCUGUUGUUCACCCACUCGUCGAACCGGCUGAAAGGGCACGAUUCUCGGUUAUAGUCCGGCCAAUUAUGAGGUAUUACAGUUUACACGGCCAUGCCUUUCCGCGUCGUGCUACGCCCACCCACUGCGAGCCGGCCUCCUUUCGAAAGCAAAGUGAUGCGUCGAAAAAGAUGAGUGUCUGCUACGGACUUAAUCAGACAGUUCCGCCUUUGAUCUGUCUCCAGGAGAAUUUUAGGAAUCUAACGACGACCUCGGACCGGGCGAAGUAUUCUUACGUGGCCACUUUUGAGCGAAUUGAGAGCUCGACAAGAUGCGUGAUAGCAAUCCACCCUACCAGGCUUUUGCGAGUGGCUAUAUGUGGACCGUUAAAACCCCUCGGGCCCCUUCAAACACUCUAUCCAUCAUUCUGCGUUCAUUCUUACAACCCAGUGACAAAUCGGGGUCAGAAGCAGUCGAGCGGGCAGGAGCCCAAACUACGUGACGUCACUACGACGCAGCCCCAACAUGGAGCUUUCGAGACCUCUGAAUCCUACCGUGUGAGUUUACAACACUCACACCCAACGAAGCGUCCCGUUCGGGGGGCUCCAUGGCUGGAAUUUCUGACGAACACCAGACUCCGUUUUAUCCUCACUAAAAAUGUACCGUGCGCCACUGCCCAACGUUCGGACGUUAGACCAGGUUCUAUACCAGCUGUCUUCAUACCAUCGGUUGAUGUAGAUGUCCGCACCAUUCUUCGCCCCAAGCCACUUACGGUGGAGUGGCCGUCAUCUCAGUGUCCUGAAGUUUAUAUGCUACUUUUUGUUGGAUACAUACGUCGCGAUAAAUUAUUUCGUAGUUCAACGUACGAACAUUACGCCUGUCACCCUAGGGAAGAUUAUUCCCCGAGCAGUUCCAGAGUGCAUGACAUGUGGAUGUUUGAUGCGGCAUGGCUGCUGGUGAGUCCUCCUUCUGUUACUAACGUAUCAUAUCAGAGUCAUAAUCAUCUGCUGCCAGCAGUCCAGAUAGGAAGUUCUCAAAGCCGCCUUGAAUGUUCAGGACCGGCAUAUUGCAUCCGCACGCUGGGGAUGCCACGUGAUACCGCAAAACACCAACACCUGUCGCGGUACGGGUCGUGGUUUAUCCGUAAUAUUCACGAUUUGCGGAAGUCACAUCUUGCAGAUUGCUACUCUGGGCACCUGUGCAUAGGCUAUCCUCUUGUAUGCACCAUGCGCCUUGUACUCGGACAAAGGGAGCUAUUACCAACUUGCCCCCCGUUGCCUGGAAAGAAGUUGCUGGAGUCGGAUGAUUUGCAUCUGCAUCCCCCUCUGGAAGUGUUGUUCCCGUUACUCUUGGUGCUCAAAGCUACGCACCUCCAGGGAGGCAGUCAUUCGGCGGUGCUGUGUUAUCAAGUGUUGUCAGCCAUACACGCUUUCCCUGUCGCUUGGCAUGGUGUCAGUGCAUACGACCUCCGACAAAUCAUAGAGAUGUUUGUUUGGUUGGGACGUAGCAAGCACAAGGGAGGAAACGGCACCUCCCUAUUUGGCGUACUUCGCUCGUGCUCUUCACGACGAAGUAGAAGGGAACGGGUACCGGUGAAAUCGCGUGGGGAGGCUAACGUCAACCAACAAAGAAGACGCCGUCAGUUAGCUACGAGAUGCCGCUAUCCAAAUCCACUACAAGUCCGUGUACACCUAGUCCCUGACCUUGGAAAAAACGAUCCUGGGGGCUCAAGCCCACGUCUAUUGACACUGGGUCGAGGUUUAGAGUGGAAGGCUGAGGCAUGCCAGUGGAAACCCUUAUUUUGCCUGGUAGGGCUGCGGUCUACGGGCCUCAUGGAUUCCAAAUUAAUGUACAAUAAAGGUACAAUACUUCAAAGAAAGAUAUUGAGACCUCCUACCCUAGAAGCGGUGAGGCUGAAAUUGGCCCAUACCCGUCUGGACAAGGAUUCCCGCUCUCCGAGGCAGAAGAGAGUAAACACCACCCUCGGUUCGGCUCCGGCAUUGUUAGCUGGAAGCGCUUGGAUUGCUAGCGGCCGGCGGGAUACGACUAGGCCGAAGAACGGAUUCUUGUCUUGCGAGUUAAUAAUAAGGGCCGCUCAAGAGACUGCGACGAACAGCCACAAUCAUGACGCGGGCAUAUCCCCAGCAAGCGCCACAGCAUCUUGCUGGGUUCGGCCUUUCCGGUACACCACAGAAGUCAGAGAGAAUCCCCCAAUAACCCGACCCCUAAGGUCUAGUGUGAUAUACCUGCCUUCGCUCCAUCUGGGAGGCUUCGAUUGCUUGACACCGACACUCCCGAAUGGAAAACUUUCCGCCAUCUACGGUGUUGUGAACACCUCGGGGCGCAAUCGCAACCUUCUUUUCCGGGCUGCCCUCGGGUGCUCGGAGACGUUCGCAAAUGUUGCAGACGCUUACGGGUCAACUCAUGCUCGCCCGGACUCCCUAUACGCUACACCGGUGGCCGGGAGAAGCAUGGCGGUUACGUUAUCAAACUCAACUACGGGGAGUACAAACGCACCCCGCGUGAAAUUCCUAGACCGCAAGGGACCCCUCAGUACCGCGGUCAGCCAGAAGUUCUUGGAAGUCCCGGUAGAUCGGCGGACUGGUUUACAUAUAUCUACCAUCACCAGCAUGAACGCAAGGCGCAUCCAGUGUGCGGGUAACUAUGGCAGCGCUUUAUGUUGCCUUCCACUCGGCCUGGUCCUGCCGACUGCACUACGCGAAGGCAGGAUUAUGUCACUGUGCCACUGGCCGAAGAGCGCGAGCGGGGCGACGGUAUUGUGUACCAAUGUACCUUCCUGCGAGCGUCCGACCGGGUGUUUUGCUACAUUCGGCCACGUCGCUUGGUGUUUGACGGUAGAUGGGCUACCUCAACGAAUCACUACGGCGUGCGAUAGAGAGUGGUCUUGCCCGUUAAGCCAGGUAAUUUUAAGAGCUGGUGUUUUAUACAUCGCGCAAACCAUAACGUCCUUACCGGAACCUAACACGACAGCGCAAGAGUCCCGCGUCAGAAACGUGCGCAUCUAUAAGCGUACCGUACCUACUUUGGCACGCAUACCACCACGGAACGUUAAAAAUGAUUUCUGCAGCACUUCUCGACGAUUAUCGGCCUACCUUGGCUUCACGAGUUCUAAUCCUUGGCCUCUAGAUCCUAAGGGAUUCAUGGUGGCAGACCUUUAUCCGCAGCUCACAGUUGGUGUAUCUGCUAUCCACUCUCAACAGCUAGCAAGUCCAUUUCUAGAUUGCCGGUCCCUAAGGAUAUUCACGAAUUCUUGGCGCCAGAAGCUCAUAUCCAUGAAAGAGAUGACUAGUAAGUCAUCAAACAGGUACUUACCCAUUUGUUCAGAGAUUCUUGACGUUCUUCGUGAGCUGUCAUAUGUGAACCGGCUUUCGCUGGGGCUACAGUUUACAGUCCCGGCGAGACAAAGGAUAACUGUGCCGGCAUCUGCUGCCCGGUUACUCAGAAUGCGGAGUCAGUUACGUAAUCUCUUAUCUAAUAAGGCACGUGUAGCGCACGUACACCAAGACCAUCAGACUGUAAAGCCAAUACGCGGAUUCGCUCUGGCCAUUCGCCAGGAUAUCCUGCGUAGCCUUGAGCGCCCUGGCCAGGACGAGCAGGAUUCCCACUGUACUCGUCCCUAUGACCACGUGGUAAGACGAGAGGCGAGUUCGUGCGGUUUCCACCCGUGCCGGCAUCUUAACGGGUCAAAGACCGCAUCACCACUAACUUCAGUCGCGAGAAAUGUCGCACCCGUCCUACACAAGACACCGAGCUUACUUUCGGCGUAUGGGCGAUGUCUAAUACGGACGCGCCCUAUCACUUAA";
    
    //traduzir_rna_proteina_single(proteina);
    //traduzir_rna_proteina(proteina);
    
}

//MRIPGAVTFVRRYIGNMNRSGRHLLSVLYVCSGGKRIPPRYMHLMELGGTPEVRTEVNRGRGPCRHTSNPTVNISRGQRLSDSKHCRGFYMPSLLSVNTHGSASDNRPAGISLWQLELGTWFTQYLFEIKLHKNVIALLSTSRLWKTRNGRLSGNYPSHPRHNCSPTSGPASNLVALRLSKAGLNPSIVVFCHRKELKLALIYVVTNLSRTGNWQIILARPWHGPALYRTKGPTTVPSFHFCSPVMDRGVLGRTKDSQDRQIAASCLPAPTVMVTRSGHTALTSGWNIDVMYKKLAGLKHLVASAFPHSPSSLFSDKIQRRVVTSGFYAPVWGFALVQDTEHRYCSTARTTSGLCDICQHSYESETDSRPCSRTRQSVQGLGSYGQRCGRPTARQTLPPEGTEPGPVGDLESGCMEFNVRLLVLGRHTLSVWCGLTRVANCSDGTSVDRLQAGELLISDVGHYVNFCLTYLPETVHPCARVGCLPTFNRIIRTFVRLVTPHLMLFVLLTAALMTQAVCWDKRYSTKRTLESLSGVHYSHNLQTTVRACGIRWMLCGFRTKFQRSANSRTCVSLMVIGRIQGKKLQATARDKMFSVEAHFQVLKVTNLSAAKLELISLGKQEAEPQQTLTLHPHIGAQYRVFSDPGVRGVPRFRQHNHYCVSNWEPGTFKETWEFAQEQNNWSFSKDEVACKIPEYRWCRIDSLRSLLYAHVYAHKLVQSTRFERLAGVGCLHLIALGPVWLTKAQTLLTLGNRVSDLKGCYIGTNLVLLHHVICSEYGCKSITVQVRCFINIATRICGVYVWSFLTTPYGTPSRQDISRVYVVPNTMPFIGSGFVRGKKALPPNCKGRSGKATLHTCEPCWLHDPILVYIVGAMNWRRVQKSVAQFSLRVSGIPGNSLRFSNPTCGGPRIYTPQSLDPACCCGREHEMDVCLVQIEIYGESKSILPKDPAGNELEHCGFPNFSCIRLFLWNHMWRVACQSPVIAYFKTVTTGAALRRNTSLSCTRLRPDAGCHLAREYIVCIGVKLPDEPLVEKSLTYTVPDVARVKLIATKVLVRECFTLHSPKGYVLDDRIDFVTDATNQETVSYMSILTCSEDVLIKRMSRVPSVCWLATFKILAPRTWRAIIRSRNAVSLLFRGGILEPRLYLCTCAYEGAVMESVHMIGHKLVGGKGAFCKTIAAKARYKECVPPIQQANSARYQNYRELPKVAKCSRVALQFGLWCEFRSLTAPPMLYARANLLSIDNFLRAPCFSSTIHDSVRLFHLQMVIKKSRRIAVCIRPPLSERRGSPCINCEGVRTASRKIITLHAKILTATIAVPGYCRADSGSPGENEYSETVQAHLAQVGATFSKGLPGTRTACWVPSPKVNLTSIPKVQYISFRRKSWELHCASQPLAVWSVDSVHWGLQRYMLPVNKKSRKVLEVNNSSIITAPNDQHSYHAIGHWTTPVVHPLVEPAERARFSVIVRPIMRYYSLHGHAFPRRATPTHCEPASFRKQSDASKKMSVCYGLNQTVPPLICLQENFRNLTTTSDRAKYSYVATFERIESSTRCVIAIHPTRLLRVAICGPLKPLGPLQTLYPSFCVHSYNPVTNRGQKQSSGQEPKLRDVTTTQPQHGAFETSESYRVSLQHSHPTKRPVRGAPWLEFLTNTRLRFILTKNVPCATAQRSDVRPGSIPAVFIPSVDVDVRTILRPKPLTVEWPSSQCPEVYMLLFVGYIRRDKLFRSSTYEHYACHPREDYSPSSSRVHDMWMFDAAWLLVSPPSVTNVSYQSHNHLLPAVQIGSSQSRLECSGPAYCIRTLGMPRDTAKHQHLSRYGSWFIRNIHDLRKSHLADCYSGHLCIGYPLVCTMRLVLGQRELLPTCPPLPGKKLLESDDLHLHPPLEVLFPLLLVLKATHLQGGSHSAVLCYQVLSAIHAFPVAWHGVSAYDLRQIIEMFVWLGRSKHKGGNGTSLFGVLRSCSSRRSRRERVPVKSRGEANVNQQRRRRQLATRCRYPNPLQVRVHLVPDLGKNDPGGSSPRLLTLGRGLEWKAEACQWKPLFCLVGLRSTGLMDSKLMYNKGTILQRKILRPPTLEAVRLKLAHTRLDKDSRSPRQKRVNTTLGSAPALLAGSAWIASGRRDTTRPKNGFLSCELIIRAAQETATNSHNHDAGISPASATASCWVRPFRYTTEVRENPPITRPLRSSVIYLPSLHLGGFDCLTPTLPNGKLSAIYGVVNTSGRNRNLLFRAALGCSETFANVADAYGSTHARPDSLYATPVAGRSMAVTLSNSTTGSTNAPRVKFLDRKGPLSTAVSQKFLEVPVDRRTGLHISTITSMNARRIQCAGNYGSALCCLPLGLVLPTALREGRIMSLCHWPKSASGATVLCTNVPSCERPTGCFATFGHVAWCLTVDGLPQRITTACDREWSCPLSQVILRAGVLYIAQTITSLPEPNTTAQESRVRNVRIYKRTVPTLARIPPRNVKNDFCSTSRRLSAYLGFTSSNPWPLDPKGFMVADLYPQLTVGVSAIHSQQLASPFLDCRSLRIFTNSWRQKLISMKEMTSKSSNRYLPICSEILDVLRELSYVNRLSLGLQFTVPARQRITVPASAARLLRMRSQLRNLLSNKARVAHVHQDHQTVKPIRGFALAIRQDILRSLERPGQDEQDSHCTRPYDHVVRREASSCGFHPCRHLNGSKTASPLTSVARNVAPVLHKTPSLLSAYGRCLIRTRPIT

//gcc -fopenmp dna.cpp -o dna.o -lstdc++; ./dna.o

