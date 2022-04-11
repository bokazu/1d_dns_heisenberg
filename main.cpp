#include "all.h"

using namespace std;

int main(){
    cout << "START 1d Heisenberg Model Calculation\n\n";
    cout << "1. Calculation Hamiltonian Matrix's Components\n\n";

    /*Hamiltonian*/
    int tot_sitenum,mat_dim,bond_num,precision;
    double *H = new double[mat_dim*mat_dim];
    vec_init(mat_dim*mat_dim,H);
    /***********Calculate Hamiltonian Matrix's Components*********/
    /***********INPUT DATA*********************/
    string M_H_settingfile_name =
        "../model_set/settingfile.txt";  //系のsite数、output用ファイル名の情報をこのfileに書いておく
    string M_H_JsetFile_name;
    string M_H_OutputFile_name;  // output用file

    ifstream M_H_settingfile(M_H_settingfile_name);

    string ftmp, Boundary_Condition;
    stringstream ss;

    make_hamiltonian(H);
    
    delete[] H;
}