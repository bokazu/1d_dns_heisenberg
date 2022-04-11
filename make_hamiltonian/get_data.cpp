#include "../all.h"

using namespace std;

void get_data(int tot_sitenum,int bond_num,std::string M_H_OutputFile_name){
    string M_H_settingfile_name =
        "../model_set/settingfile.txt";  //系のsite数、output用ファイル名の情報をこのfileに書いておく
    string M_H_jsetFile_name;

    ifstream M_H_Settingfile(M_H_settingfile_name);

    string ftmp, Boundary_Condition;
    stringstream ss;
    int fcount = 0;
}