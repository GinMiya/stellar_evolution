#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <sys/stat.h>
#include <omp.h>

using namespace std;

// 実行するプログラムの選択どれか1つ

#define DENSITY_TEMPERATURE
// #define MASS_TMAX
// #define MINIMUM_MASS

#define MAKE_PICTURE

#define G     (6.67*pow(10.0, -8.0))     // dyn cm^2 g^-2
#define kmp   (8.31433*pow(10.0, 7.0))   // erg g^-1 K^-1
#define phi   (16.15)
#define h     (6.6261*pow(10.0, -27.0))  // erg s
#define me    (9.1094*pow(10.0, -28.0))  // g
#define mp    (1.6726*pow(10.0, -24.0))  // g
#define M_sun (1.989*pow(10.0, 33.0))    // g

#ifdef DENSITY_TEMPERATURE
  #define min_rho 2    // 密度の最小order
  #define max_rho 7    // 密度の最大order
  #define bin   1000   // 1order当たりのサンプリング数

  void initialize_T(double T[], double rho[]){
    // 密度サンプリングは10^2~10^7までとって各オーダーで100サンプリングする
    int index;
    for(int i=0; i<=(max_rho-min_rho)*bin; i++){
      T[i] = 0.0;
    }
    for(int i=min_rho; i<max_rho; i++){
      for(int j=0; j<bin; j++){
        index = (i-min_rho)*bin+j;
        rho[index] = pow(10.0, i) + pow(10.0, i+1)*j/bin;
      }
    }
  }
  double calc_mu_H_He_Z(double X, double Y, double Z){
    // // 電離ガスの平均分子量の計算
    double mu_H = 2.0*X;
    double mu_He = 3.0*Y/4.0;
    double mu_Z = 0.5*Z;
    return 1.0/(mu_H+mu_He+mu_Z);
  }
  double calc_mu_electron(double X, double Y, double Z){
    // // 電離ガスの平均分子量の計算
    double mu_H_e = X;
    double mu_He_e = 0.5*Y;
    double mu_Z_e = 0.5*Z;
    return 1.0/(mu_H_e+mu_He_e+mu_Z_e);
  }
  double calc_ne(double rho, double mu_e){
    return rho/(mp*mu_e);
  }
  double calc_degenerate_P(double ne){
    double factor = (h*h)*pow(3.0/M_PI, (2.0/3.0))/(20.0*me);
    return factor * pow(ne, (5.0/3.0));
  }
  double calc_T(double rho, double M, double mu, double ne){
    double numer, denom;
    numer = 4.0*M_PI*(G*G*G)*(M*M)*(rho*rho*rho*rho)/(phi*phi);
    numer = pow(numer, (1.0/3.0));
    denom = kmp*rho/(mu);
    numer -= calc_degenerate_P(ne);
    return numer/denom;
  }
  void set_T_array(double T[], double rho_c[], double M, double mu, double mu_e, int Nbin){
    double ne;
    #pragma omp parallel for
    for(int i=0; i<=Nbin; i++){
      ne = calc_ne(rho_c[i], mu_e);
      T[i] = calc_T(rho_c[i], M, mu, ne);
    }
    return;
  }
  void outputfile_1(double T_c[], double rho_c[], int Nbin, int id){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_HHeZ");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/dens_temp_%d.dat", id);
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<=Nbin; i++){
      outputfile << rho_c[i] << "\t" << T_c[i] << endl;
    }
    outputfile .close();
  }

  double calc_mu_C_O(double X_C, double X_O){
    double mu_C = 0.5+X_C/12.0;
    double mu_O = 0.5+X_O/16.0;
    return 1.0/(mu_C+mu_O);
  }
  double calc_mu_e_C_O(double X_C, double X_O){
    double mu_C = X_C*0.5;
    double mu_O = X_O*0.5;
    return 1.0/(mu_C+mu_O);
  }
  void outputfile_2(double T_c[], double rho_c[], int Nbin, int id){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_CO");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/dens_temp_%d.dat", id);
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<=Nbin; i++){
      outputfile << rho_c[i] << "\t" << T_c[i] << endl;
    }
    outputfile .close();
  }

  int main(void){
    // ================================================= //
    // H,He,Zの完全電離ガスにおける中心密度-温度関係の作成   //
    // ================================================= //
    cout << "1. density-temperature relation of H,He,Z star." << endl;
    double X, Y, Z; // H,He,metalの存在比
    cout << "Input X, Y, Z: ";
    cin >> X >> Y >> Z;
    double mu = calc_mu_H_He_Z(X, Y, Z);
    double mu_e = calc_mu_electron(X, Y, Z);
    double M_ch = 1.46*pow(2.0/mu_e, 2.0)*M_sun;
    cout << "Average molecular wight: " << mu << endl;
    cout << "Average electron wight: " << mu_e << endl;
    cout << "Chandrasekhar mass: " << M_ch/M_sun << endl;
    int N_bin = (max_rho-min_rho)*bin;
    double T_c[N_bin+1], rho_c[N_bin+1];
    initialize_T(T_c, rho_c);

    // 質量変えて探査
    double M;
    for(int i=1; i<=45; i++){
      M = 0.1*i*M_sun;
      set_T_array(T_c, rho_c, M, mu, mu_e, N_bin);
      outputfile_1(T_c, rho_c, N_bin, i);
    }


    // ================================================= //
    //    C,Oの完全電離ガスにおける中心密度-温度関係の作成   //
    // ================================================= //
    cout << "2.density-temperature relation of C,O star." << endl;
    double X_C, X_O; // C,Oの存在比
    cout << "Input X_C, X_O: ";
    cin >> X_C >> X_O;
    mu = calc_mu_C_O(X_C, X_O);
    mu_e = calc_mu_e_C_O(X_C, X_O);
    M_ch = 1.46*pow(2.0/mu_e, 2.0)*M_sun;
    cout << "Average molecular wight: " << mu << endl;
    cout << "Average electron wight: " << mu_e << endl;
    cout << "Chandrasekhar mass: " << M_ch/M_sun << endl;
    initialize_T(T_c, rho_c);
    M = 0.0;
    for(int i=1; i<=20; i++){
      M = 0.1*i*M_sun;
      set_T_array(T_c, rho_c, M, mu, mu_e, N_bin);
      outputfile_2(T_c, rho_c, N_bin, i);
    }

    #ifdef MAKE_PICTURE
      FILE *fpg;
      fpg = popen("gnuplot", "w");
      if(fpg == NULL) return -1;
        fputs("load 'gnuplot/pic.gp'\n",fpg);
        fputs("load 'gnuplot/pic_CO.gp'\n",fpg);
        fputs("load 'gnuplot/pic_all.gp'\n",fpg);
        fputs("load 'gnuplot/pic_all_CO.gp'\n",fpg);
        fflush(fpg);
      pclose(fpg);
    #endif
    return 0;
  }
#endif

#ifdef MASS_TMAX
  #define min_M 0.1    // 最小質量
  #define max_M 10     // 最大質量
  #define bin   100   // 1order当たりのサンプリング数

  void initialize_T_M(double T[], double M[]){
    for(int i=0; i<=(max_M-min_M)*bin; i++){
      T[i] = 0.0;
      M[i] = M_sun*min_M + M_sun*((max_M-min_M))*i/((max_M-min_M)*bin);
    }
  }

  double calc_mu_H_He_Z(double X, double Y, double Z){
    // // 電離ガスの平均分子量の計算
    double mu_H = 2.0*X;
    double mu_He = 3.0*Y/4.0;
    double mu_Z = 0.5*Z;
    return 1.0/(mu_H+mu_He+mu_Z);
  }
  double calc_mu_electron(double X, double Y, double Z){
    // // 電離ガスの平均分子量の計算
    double mu_H_e = X;
    double mu_He_e = 0.5*Y;
    double mu_Z_e = 0.5*Z;
    return 1.0/(mu_H_e+mu_He_e+mu_Z_e);
  }
  double calc_T_max_independent(double M, double mu, double mu_e){
    double numer, denom;
    numer = 4.0*M_PI*(G*G*G)*(M*M)/(phi*phi);
    numer = pow(numer, (2.0/3.0));
    denom = 4.0*(kmp/mu)*(h*h)/(20.0*me)*pow(3.0/M_PI, (2.0/3.0))/pow(mp*mu_e, (5.0/3.0));
    return numer/denom;
  }
  void set_T_max_array(double T[], double M[], double mu, double mu_e, int Nbin){
    double ne, rho_c;
    #pragma omp parallel for
    for(int i=0; i<=Nbin; i++){
      T[i] = calc_T_max_independent(M[i], mu, mu_e);
    }
    return;
  }
  void outputfile_1(double T_c[], double M[], int Nbin){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_tmax_HHeZ");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/tmax.dat");
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<=Nbin; i++){
      outputfile << M[i] << "\t" << T_c[i] << endl;
    }
    outputfile .close();
  }

  double calc_mu_C_O(double X_C, double X_O){
    double mu_C = 0.5+X_C/12.0;
    double mu_O = 0.5+X_O/16.0;
    return 1.0/(mu_C+mu_O);
  }
  double calc_mu_e_C_O(double X_C, double X_O){
    double mu_C = X_C*0.5;
    double mu_O = X_O*0.5;
    return 1.0/(mu_C+mu_O);
  }
  void outputfile_2(double T_c[], double M[], int Nbin){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_tmax_CO");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/tmax.dat");
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<=Nbin; i++){
      outputfile << M[i] << "\t" << T_c[i] << endl;
    }
    outputfile .close();
  }

  int main(void){
    // ================================================= //
    // H,He,Zの完全電離ガスにおける質量-最高温度関係の作成   //
    // ================================================= //
    cout << "3-1.Max temperature-Mass relation of H,He,Z star." << endl;
    double X, Y, Z;
    cout << "Input X, Y, Z: ";
    cin >> X >> Y >> Z;
    double mu = calc_mu_H_He_Z(X, Y, Z);
    double mu_e = calc_mu_electron(X, Y, Z);
    double M_ch = 1.46*pow(2.0/mu_e, 2.0)*M_sun;
    cout << "Average molecular wight: " << mu << endl;
    cout << "Average electron wight: " << mu_e << endl;
    cout << "Chandrasekhar mass: " << M_ch/M_sun << endl;

    int N_bin = (max_M-min_M)*bin;
    double T_c[N_bin+1], M[N_bin+1];
    initialize_T_M(T_c, M);
    set_T_max_array(T_c, M, mu, mu_e, N_bin);
    outputfile_1(T_c, M, N_bin);


    // ================================================= //
    //    C,Oの完全電離ガスにおける質量-最高温度関係の作成   //
    // ================================================= //
    cout << "3-2.Max temperature-Mass relation of C,O star." << endl;
    double X_C, X_O;
    cout << "Input X_C, X_O: ";
    cin >> X_C >> X_O;
    mu = calc_mu_C_O(X_C, X_O);
    mu_e = calc_mu_e_C_O(X_C, X_O);
    M_ch = 1.46*pow(2.0/mu_e, 2.0)*M_sun;
    cout << "Average molecular wight: " << mu << endl;
    cout << "Average electron wight: " << mu_e << endl;
    cout << "Chandrasekhar mass: " << M_ch/M_sun << endl;

    initialize_T_M(T_c, M);
    set_T_max_array(T_c, M, mu, mu_e, N_bin);
    outputfile_2(T_c, M, N_bin);

    #ifdef MAKE_PICTURE
    FILE *fpg;
    fpg = popen("gnuplot", "w");
    if(fpg == NULL) return -1;
      fputs("load 'gnuplot/pic_tmax.gp'\n",fpg);
      fputs("load 'gnuplot/pic_tmax_CO.gp'\n",fpg);
      fflush(fpg);
    pclose(fpg);
    #endif

    return 0;
  }
#endif

#ifdef MINIMUM_MASS
  double calc_mu_H_He_Z(double X, double Y, double Z){
    double mu_H = 2.0*X;
    double mu_He = 3.0*Y/4.0;
    double mu_Z = 0.5*Z;
    return 1.0/(mu_H+mu_He+mu_Z);
  }
  double calc_mu_electron(double X, double Y, double Z){
    double mu_H_e = X;
    double mu_He_e = 0.5*Y;
    double mu_Z_e = 0.5*Z;
    return 1.0/(mu_H_e+mu_He_e+mu_Z_e);
  }
  double calc_mu_metal(int A){
    double mu = 0.5+1.0/A;
    return 1.0/mu;
  }
  void calc_mu_T_Mch(int i, double mu[], double mu_e[], double T[], double M_ch[]){
    switch(i){
      case 0:
        cout << "== Composition: H, He, Z ==" << endl;
        double X, Y, Z;
        cout << "Input X, Y, Z: ";
        cin >> X >> Y >> Z;
        mu[0] = calc_mu_H_He_Z(X, Y, Z);
        mu_e[0] = calc_mu_electron(X, Y, Z);
        T[0] = 2.0*pow(10.0, 7.0);
        M_ch[0] = 1.46*pow(2.0/mu_e[0], 2.0)*M_sun;
        cout << "\tMu: " << mu[0] << " Mu_e: " << mu_e[0] << endl;
        break;
      case 1:
        cout << "== Composition: He ==" << endl;
        mu[1] = calc_mu_metal(4.0);
        mu_e[1] = 2.0;
        T[1] = 2.0*pow(10.0, 8.0);
        M_ch[1] = 1.46*pow(2.0/mu_e[1], 2.0)*M_sun;
        cout << "\tMu: " << mu[1] << " Mu_e: " << mu_e[1] << endl;
        break;
      case 2:
        cout << "== Composition: C ==" << endl;
        mu[2] = calc_mu_metal(12.0);
        mu_e[2] = 2.0;
        T[2] = 7.0*pow(10.0, 8.0);
        M_ch[2] = 1.46*pow(2.0/mu_e[2], 2.0)*M_sun;
        cout << "\tMu: " << mu[2] << " Mu_e: " << mu_e[2] << endl;
        break;
      case 3:
        cout << "== Composition: Ne ==" << endl;
        mu[3] = calc_mu_metal(20.2);
        mu_e[3] = 2.0;
        T[3] = pow(10.0, 9.0);
        M_ch[3] = 1.46*pow(2.0/mu_e[3], 2.0)*M_sun;
        cout << "\tMu: " << mu[3] << " Mu_e: " << mu_e[3] << endl;
        break;
      case 4:
        cout << "== Composition: O ==" << endl;
        mu[4] = calc_mu_metal(16.0);
        mu_e[4] = 2.0;
        T[4] = 2.0*pow(10.0, 9.0);
        M_ch[4] = 1.46*pow(2.0/mu_e[4], 2.0)*M_sun;
        cout << "\tMu: " << mu[4] << " Mu_e: " << mu_e[4] << endl;
        break;
      default:
        break;
    }
  }
  void set_mu_mue_T(double mu[], double mu_e[], double T[], double M_ch[]){
    for(int i=0; i<5; i++) calc_mu_T_Mch(i, mu, mu_e, T, M_ch);
  }
  double calc_min_mass(double T, double mu, double mu_e){
    double numer, denom, M_min;
    denom = 4.0*M_PI*(G*G*G)/(phi*phi);
    denom = pow(denom, (2.0/3.0));
    numer = 4.0*(kmp/mu)*(h*h)/(20.0*me)*pow(3.0/M_PI, (2.0/3.0))/pow(mp*mu_e, (5.0/3.0));
    M_min = pow(numer*T/denom, (3.0/4.0));
    return M_min;
  }
  void set_min_mass(double T[], double M[], double mu[], double mu_e[]){
    for(int i=0; i<5; i++){
      M[i] = calc_min_mass(T[i], mu[i], mu_e[i]);
    }
  }
  void outputfile(double T[], double M[]){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_mmin");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/mmin_t.dat");
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<5; i++){
      outputfile << M[i] << "\t" << T[i] << endl;
    }
    outputfile .close();
  }

  #define min_M 0.1
  #define max_M 10.0
  #define bin   100
  void initialize_T_M(double T[], double M[]){
    for(int i=0; i<=(max_M-min_M)*bin; i++){
      T[i] = 0.0;
      M[i] = M_sun*min_M + M_sun*((max_M-min_M))*i/((max_M-min_M)*bin);
    }
  }
  double calc_T_max_independent(double M, double mu, double mu_e){
    double numer, denom;
    numer = 4.0*M_PI*(G*G*G)*(M*M)/(phi*phi);
    numer = pow(numer, (2.0/3.0));
    denom = 4.0*(kmp/mu)*(h*h)/(20.0*me)*pow(3.0/M_PI, (2.0/3.0))/pow(mp*mu_e, (5.0/3.0));
    return numer/denom;
  }
  void set_T_max_array(double T[], double M[], double mu, double mu_e, int Nbin){
    double ne, rho_c;
    #pragma omp parallel for
    for(int i=0; i<=Nbin; i++){
      T[i] = calc_T_max_independent(M[i], mu, mu_e);
    }
    return;
  }
  void outputfile_1(double T_c[], double M[], int Nbin, int id){
    //保存用フォルダの作成
    char dirname[50], _dirname[50], filename[50];
    sprintf(dirname,"data_mmin");
    strcpy(_dirname,dirname);
    mkdir(dirname,0777);
    chmod(dirname,0777);

    sprintf(filename,"/tmax_%d.dat", id);
    strcat(_dirname,filename);
    ofstream outputfile(_dirname);
    for(int i=0; i<=Nbin; i++){
      if(T_c[i] == 0.0) continue;
      outputfile << M[i] << "\t" << T_c[i] << endl;
    }
    outputfile .close();
  }

  int main(void){
    cout << "4.Minumum Mass Burning temperature relation calculator." << endl;
    double mu[5], mu_e[5], T[5], M_min[5], M_ch[5];
    set_mu_mue_T(mu, mu_e, T, M_ch);
    set_min_mass(T, M_min, mu, mu_e);
    outputfile(T, M_min);

    cout << "4-0.Max temperature-mass relation of some Composition." << endl;
    int N_bin = (max_M-min_M)*bin;
    double T_c[N_bin+1], M[N_bin+1];
    for(int i=0; i<5; i++){
      initialize_T_M(T_c, M);
      set_T_max_array(T_c, M, mu[i], mu_e[i], N_bin);
      outputfile_1(T_c, M, N_bin, i);
    }
    #ifdef MAKE_PICTURE
    FILE *fpg;
    fpg = popen("gnuplot", "w");
    if(fpg == NULL) return -1;
      fputs("load 'gnuplot/pic_tmax_minm.gp'\n",fpg);
      fflush(fpg);
    pclose(fpg);
    #endif

    return 0;
  }
#endif
