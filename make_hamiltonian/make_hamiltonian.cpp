/*M_H = Make_Hamiltonian*/
#include "../all.h"

using namespace std;

void make_hamiltonian(int mat_dim, int tot_site_num,
                      std::string M_H_JsetFile_name,
                      std::string M_H_OutputFile_name, int precision,
                      std::string Boundary_Condition, double *H)
{
    int bond_num;

    /*jset.txtからのbondごとの相互作用情報の取得*/
    /*bond数の取得*/
    ifstream M_H_JsetFile(M_H_JsetFile_name);
    if (!(M_H_JsetFile))
    {
        cerr << "Could not open the file(line 10) - '" << M_H_JsetFile_name
             << "'" << endl;
    }
    if (Boundary_Condition == "y")
    {
        bond_num = tot_site_num;
    }
    else
    {
        bond_num = tot_site_num - 1;
    }

    double *J = new double[bond_num];
    std::cout << "i"
              << "  "
              << "i+1"
              << ":  "
              << " J[i]      " << endl;
    for (int i = 0; i < bond_num; i++)
    {
        J[i] = 0.;
        M_H_JsetFile >> J[i];
        std::cout << i << "   " << i + 1 << "  :  " << J[i] << endl;
    }

    M_H_JsetFile.close();

    if (Boundary_Condition == "y")
    {
        for (int site_num = 0; site_num < tot_site_num; site_num++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                spm(j, site_num, tot_site_num, mat_dim, H, J);
                smp(j, site_num, tot_site_num, mat_dim, H, J);
                szz(j, site_num, tot_site_num, mat_dim, H, J);
            }
        }
    }
    else if (Boundary_Condition == "n")
    {
        for (int site_num = 0; site_num < tot_site_num - 1; site_num++)
        {
            for (int j = 0; j < mat_dim; j++)
            {
                spm(j, site_num, tot_site_num, mat_dim, H, J);
                smp(j, site_num, tot_site_num, mat_dim, H, J);
                szz(j, site_num, tot_site_num, mat_dim, H, J);
            }
        }
    }
    else
    {
        cout << "ERROR : Maybe inputed other than \"y\" and \"n\" " << endl;
    }

    // /*OUTPUT HAMILTONIAN*/
    ofstream M_H_Output(M_H_OutputFile_name);

    printmat(mat_dim, precision, H);
    fprintmat(M_H_Output, mat_dim, precision, H);

    M_H_Output.close();
    delete[] J;
}