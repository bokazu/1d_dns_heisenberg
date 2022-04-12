#include "all.h"

using namespace std;

int main()
{
    cout << "START 1d Heisenberg Model Calculation\n\n";
    cout << "1. Calculation Hamiltonian Matrix's Components\n\n";

    /*Hamiltonian*/
    int tot_site_num, precision;

    /***********************************Calculate Hamiltonian Matrix's
     * Components************************************************/

    /***********************************************************INPUT
     * DATA*******************************************************/
    string M_H_settingfile_name = "./model_set/settingfile.txt";
    ifstream if_M_H_Settingfile(M_H_settingfile_name);
    if (!(if_M_H_Settingfile))
    {
        cerr << "Could not open the file" << M_H_settingfile_name << "'"
             << endl;
    }

    string M_H_OutputFile_name, M_H_JsetFile_name, Boundary_Condition,
        D_L_OutputFile_name;
    get_data(if_M_H_Settingfile, tot_site_num, M_H_OutputFile_name,
             M_H_JsetFile_name, D_L_OutputFile_name, Boundary_Condition,
             precision);

    std::cout << "/************************************************************"
                 "***************************"
              << "INPUT DATA"
              << "*************************************************************"
                 "**************************/\n";
    std::cout << "tot_site_num        = " << tot_site_num << endl;
    std::cout << "M_H_OutputFile_name = " << M_H_OutputFile_name << endl;
    std::cout << "M_H_JsetFile_name   = " << M_H_JsetFile_name << endl;
    std::cout << "D_L_OutputFile_name = " << D_L_OutputFile_name << endl;
    std::cout << "Boundary Condition  = " << Boundary_Condition << endl;
    std::cout << "precision           = " << precision << endl;
    std::cout << "/************************************************************"
                 "***************************"
              << "*************************************************************"
                 "************************************************/\n";

    /*******************************************************************************************************************************/

    /************************************Calculate Hamiltonian Matrix
     * Components****************************************************/
    int mat_dim = 1 << tot_site_num;
    double *H = new double[mat_dim * mat_dim];
    double *eigen_value = new double[mat_dim];
    vec_init(mat_dim * mat_dim, H);
    make_hamiltonian(mat_dim, tot_site_num, M_H_JsetFile_name,
                     M_H_OutputFile_name, precision, Boundary_Condition, H);
    /*******************************************************************************************************************************/

    /*********************************************************Use Lanczos
     * Method****************************************************/
    lanczos(mat_dim, H, eigen_value, D_L_OutputFile_name, precision);

    delete[] H;
    delete[] eigen_value;
}