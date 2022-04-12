#include "../all.h"

using namespace std;

lapack_int LAPACKE_dstev(int matrix_order, char jobz, lapack_int n, double *d,
                         double *e, double *z, lapack_int ldz);

void lanczos(int mat_dim, double *H, double *eigen_value,
             std::string D_L_OutputFile_name, int precision)
{
    std::cout << "/************************************************************"
                 "***************************"
              << "Lanczos Method"
              << "*************************************************************"
                 "**************************/\n";

    ofstream of_D_L_Outputfile(D_L_OutputFile_name);
    if (!(of_D_L_Outputfile))
    {
        cerr << "Could not open the file(line 10) - '" << D_L_OutputFile_name
             << "'" << endl;
    }
    // setting Initial vector & standarbilization
    double **u = new double *[mat_dim];

    double err = 1.0e-16;
    double info_sdz;
    for (int i = 0; i < mat_dim; i++)
    {
        u[i] = new double[mat_dim];
    }
    srand(time(NULL));
    for (int i = 0; i < mat_dim; i++)
    {
        for (int j = i; j < mat_dim; j++)
        {
            u[i][j] = 0.0;
            u[j][i] = 0.0;
        }
    }

    /*初期ベクトルの決定*/
    random_device rnd;
    mt19937 mt(rnd());  //メルセンヌ・ツイスタを利用
    uniform_real_distribution<> rand01(0, 1);
    for (int i = 0; i < mat_dim; i++)
    {
        u[0][i] = rand01(mt);
    }
    sdz(mat_dim, err, u[0]);

    // of_D_L_Outputfile << "u[0] = " << endl;
    // fprintvec(of_D_L_Outputfile, precision, mat_dim, u[0]);

    double *v = new double[mat_dim];
    vec_init(mat_dim, v);
    double *alpha =
        new double[mat_dim];  // Insert diagonal elements of tridiagonal matrix
    vec_init(mat_dim, alpha);
    double *beta =
        new double[mat_dim -
                   1];  // Insert subdiagonal elements of tridiagonal matrix
    vec_init(mat_dim - 1, beta);
    // Insert eigenvalue when k == even & odd
    double *eigenv_even = new double[mat_dim];
    vec_init(mat_dim, eigenv_even);
    double *eigenv_odd = new double[mat_dim];
    vec_init(mat_dim, eigenv_odd);
    // Insert eigenvector when k == even & odd
    double *eigenvec_even = new double[mat_dim];
    vec_init(mat_dim, eigenvec_even);
    double *eigenvec_odd = new double[mat_dim];
    vec_init(mat_dim, eigenvec_odd);
    // Use as lapack argument. d = alpha, e = beta
    double *d = new double[mat_dim];
    vec_init(mat_dim, d);
    double *e = new double[mat_dim - 1];
    vec_init(mat_dim - 1, e);

    double beta_pow2 = 0.;
    double eps = 1.0;
    int check = 1;  //"1" means eps > 1.0e-16. "0" means eps < 1.0e-16
    bool checker = true;
    int count = 0;
    int info;

    for (int k = 0; k < mat_dim; k++)
    {
        // cout << "count = " << k << endl;
        if (k > 0) count++;
        vec_init(mat_dim, v);
        if (checker)
        {
            if (k == mat_dim - 1)
            {
                // calculate v[i] = Au0(k)
                cblas_dgemv(CblasRowMajor, CblasNoTrans, mat_dim, mat_dim, 1.0,
                            H, mat_dim, u[k], 1, 0.0, v, 1);
                // calculate alpha & beta
                alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
            }
            else
            {
                // calculate v[i] = Au0(k)
                cblas_dgemv(CblasRowMajor, CblasNoTrans, mat_dim, mat_dim, 1.0,
                            H, mat_dim, u[k], 1, 0.0, v, 1);
                if (k == 0)
                {
                    alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
                    cblas_daxpy(mat_dim, -alpha[k], u[k], 1, v, 1);
                    beta[k] = cblas_dnrm2(mat_dim, v, 1);
                    cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
                    cblas_dcopy(mat_dim, v, 1, u[k + 1], 1);
                    sdz(mat_dim, err, u[k + 1]);
                }
                else
                {
                    alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
                    cblas_daxpy(mat_dim, -beta[k - 1], u[k - 1], 1, v, 1);
                    cblas_daxpy(mat_dim, -alpha[k], u[k], 1, v, 1);
                    beta[k] = cblas_dnrm2(mat_dim, v, 1);
                    cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
                    cblas_dcopy(mat_dim, v, 1, u[k + 1], 1);
                    sdz(mat_dim, err, u[k + 1]);
                }
            }

            // calculate eigenvalue of A(k)
            cblas_dcopy(mat_dim, alpha, 1, d, 1);
            cblas_dcopy(mat_dim - 1, beta, 1, e, 1);
            if (k % 2 == 0)
            {
                if (k == mat_dim - 1)
                {
                    info = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, d, e,
                                         eigenvec_even, mat_dim);
                    cblas_dcopy(mat_dim, d, 1, eigenv_even, 1);
                }
                else
                {
                    info = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, d, e,
                                         eigenvec_even, k + 2);
                    cblas_dcopy(mat_dim, d, 1, eigenv_even, 1);
                }
            }
            else
            {
                if (k == mat_dim - 1)
                {
                    info = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, d, e,
                                         eigenvec_odd, mat_dim);
                    cblas_dcopy(mat_dim, d, 1, eigenv_odd, 1);
                }
                else
                {
                    info = LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, d, e,
                                         eigenvec_odd, k + 2);
                    cblas_dcopy(mat_dim, d, 1, eigenv_odd, 1);
                }
            }

            // check errors each eigenvalue of groundstates
            if (k > 0)
            {
                eps = abs(eigenv_even[0] - eigenv_odd[0]);
                if (eps > err)
                {
                    checker = true;
                }
                else if (eps < err)
                {
                    checker = false;
                    // cout << "break at count = " << k << endl;
                }
            }
            // cout << "count = " << k << endl;
        }
        else
        {
            cout << "break at" << k << endl;
            break;
        }
    }
    if (count % 2 == 0)
        cblas_dcopy(mat_dim, eigenv_even, 1, eigen_value, 1);
    else
        cblas_dcopy(mat_dim, eigenv_odd, 1, eigen_value, 1);
    printf("my eigen value = \n");
    printvec(mat_dim, precision, eigen_value);
    of_D_L_Outputfile << "count = ";
    of_D_L_Outputfile << count << endl;
    of_D_L_Outputfile << endl;
    // fprintf(file, "\n");
    of_D_L_Outputfile << endl;
    of_D_L_Outputfile << "eigen value = " << endl;

    fprintvec_col(of_D_L_Outputfile, mat_dim, precision, eigen_value);

    std::cout << "end\n";
    for (int i = 0; i < mat_dim; i++)
    {
        delete[] u[i];
    }

    delete[] u;
    delete[] v;
    delete[] alpha;
    delete[] beta;
    delete[] eigenv_even;
    delete[] eigenv_odd;
    delete[] eigenvec_even;
    delete[] eigenvec_odd;
    delete[] d;
    delete[] e;
}