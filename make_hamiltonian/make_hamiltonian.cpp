/*M_H = Make_Hamiltonian*/
#include "../all.h"

using namespace std;

void make_hamiltonian(double *H)
{
    int tot_site_num, bond_num,precision;
    //<i|H|j>
    //|0> = |00>,|1> = |01>,|2> = |10>, |3> = |11> (|fermion2,fermion1>)
    string M_H_settingfile_name =
        "../model_set/settingfile.txt";  //系のsite数、output用ファイル名の情報をこのfileに書いておく
    string M_H_JsetFile_name;
    string M_H_OutputFile_name;  // output用file

    ifstream M_H_settingfile(M_H_settingfile_name);

    string ftmp, Boundary_Condition;
    stringstream ss;
    int fcount = 0;

    std::cout << "/***************************************************************************************/\n" 
    << "INPUT DATA" << "/***************************************************************************************/\n" ;
    while (!M_H_settingfile.eof())
    {
        getline(M_H_settingfile, ftmp);
        ss << ftmp;
        std::cout << "ftmp=" << ftmp << endl;
        if (fcount == 0)
        {
            ss >> tot_site_num;
            std::cout << "tot_site_num=" << tot_site_num << endl;
        }
        else if (fcount == 1)
        {
            M_H_OutputFile_name = ftmp;
            std::cout << "M_H_OutputFile_name = " << M_H_OutputFile_name << endl;
        }
        else if (fcount == 2)
        {
            M_H_JsetFile_name = ftmp;
            std::cout << "M_H_JsetFile_name : " << M_H_JsetFile_name << endl;
        }
        else if (fcount == 3)
        {
            Boundary_Condition = ftmp;
            std::cout << "Boundary Condition : " << Boundary_Condition << endl;
        }
        else if(fcount == 4)
        {
            ss >> precision;
            std::cout << "precision = " << precision << endl;
        }
        fcount += 1;
    }
    M_H_settingfile.close();
    std::cout << "/***************************************************************************************/\n";


    int dim = 1 << tot_site_num;
    std::cout << "dim=" << dim << endl;
    double *H = new double[dim * dim];
    // Hamiltonianの初期化(全要素0で埋める)
    for (int i = 0; i < dim * dim; i++) H[i] = 0.;

    /*jset.txtからのbondごとの相互作用情報の取得*/
    /*bond数の取得*/
    ifstream M_H_JsetFile(M_H_JsetFile_name);
    if (Boundary_Condition == "y")
    {
        bond_num = tot_site_num;
    }
    else
    {
        bond_num = tot_site_num - 1;
    }

    double *J = new double[bond_num];
    for(int i=0;i<bond_num;i++){
        J[i] = 0.;
        M_H_JsetFile >> J[i];
        std::cout << i << "   " << i+1 << "  :  " << "J[" << i << "]" << endl;
    }


    M_H_JsetFile.close();

    if (Boundary_Condition == "y")
    {
        for (int site_num = 0; site_num < tot_site_num; site_num++)
        {
            for (int j = 0; j < dim; j++)
            {
                spm(j, site_num, tot_site_num, dim, H,J);
                smp(j, site_num, tot_site_num, dim, H,J);
                szz(j, site_num, tot_site_num, dim, H,J);
            }
        }
    }
    else if (Boundary_Condition == "n")
    {
        for (int site_num = 0; site_num < tot_site_num - 1; site_num++)
        {
            for (int j = 0; j < dim; j++)
            {
                spm(j, site_num, tot_site_num, dim, H,J);
                smp(j, site_num, tot_site_num, dim, H,J);
                szz(j, site_num, tot_site_num, dim, H,J);
            }
        }
    }

    // /*OUTPUT HAMILTONIAN*/
    ofstream M_H_Output(M_H_OutputFile_name);

    printmat(dim,precision, H);
    fprintmat(M_H_Output, dim,precision, H);


    M_H_Output.close();
    delete[] H;
    delete[] J;
}