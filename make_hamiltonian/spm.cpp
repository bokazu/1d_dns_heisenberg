#include "../all.h"

using namespace std;
// const size_t site = 4;  // size数を変えるごとに変更が必要

void spm(int j, int site_num, int max_site_num, int dim,
         double *H,double *J)  // bra_iはsiteを変更するごとに手で変える必要あり
{
    // cout << "/************************start spm***************************/"<<endl;
    // cout << "max_site_num =" << max_site_num << endl;
    // cout << "mat_dim      = " << dim << endl; 
    boost::dynamic_bitset<> ket_j(max_site_num, j);
    boost::dynamic_bitset<> ket_j1(max_site_num, j);
    bool bit_check0, bit_check1;

    // Point A
    if (site_num != max_site_num - 1)
    {
        bit_check0 = ket_j.test(site_num);
        bit_check1 = ket_j.test(site_num + 1);
    }
    else
    {
        bit_check0 = ket_j.test(site_num);
        bit_check1 = ket_j.test(0);
    }

    if ((bit_check0 == false) && (bit_check1 == true))
    {
        // Point B
        if (site_num != max_site_num - 1)
        {
            ket_j1.flip(site_num + 1);
            ket_j1.flip(site_num);
        }
        else
        {
            ket_j1.flip(0);
            ket_j1.flip(site_num);
        }

        // Point C
        int i = ket_j1.to_ulong();
        // Hamiltonianに代入
        H[i + j * dim] += 0.5 * J[site_num];
        // cout << "H[" << j << "+" << i << "*dim] = " << H[j + i * dim] <<
        // endl;
    }
}
