// lab789_alekseeva.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <locale.h>
#include <iomanip>
using namespace std;

void BubbleSort(double* arr, int N) {
    double temp;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N - 1; j++) {
            if (arr[j] > arr[j + 1]) {
                temp = arr[j];
                arr[j] = arr[j + 1];
                arr[j + 1] = temp;
            }
        }
    }
}

void criterionShapiroWilk(double* arr, int N) {
    double standard_deviation = 0, average_value = 0;
    for (int i = 0; i < N; i++) {
        average_value += arr[i];
    }
    average_value = (1.0 / N) * average_value;

    for (int i = 0; i < N; i++) {
        standard_deviation += pow(arr[i] - average_value, 2);
    }

    int m = N / 2;
    double a0 = (0.899 / pow(N - 2.4, 0.4162)) - 0.02;
    double* z = new double[N];
    double* a = new double[N];
    double B = 0;
    for (int i = 0; i < m; i++) {
        z[i] = (N - (2.0 * (i + 1)) + 1.0) / (N - 0.5);
        a[i] = a0 * (z[i] + (1483 / pow(3 - z[i], 10.845)) + (pow(71.610, -10) / pow(1.1 - z[i], 8.26)));
        B += a[i] * (arr[N - 1 - i] - arr[i]);
    }

    B = pow(B, 2);
    double W = 0;
    W = (1 - 0.6695 / pow(N, 0.6518)) * (standard_deviation / B);
    cout << "Критерий Шапиро-Уилка" << endl;
    if (W > 1) {
        cout << "W=" << W << " нулевая гипотеза принимается" << endl;
    }
    else {
        cout << "W=" << W << " нулевая гипотеза отклоняется" << endl;
    }
    delete[] z, a;
}

void criterionLinMudholkar(double* arr, int N) {
    double standard_deviation = 0, average_value_arr = 0;
    for (int i = 0; i < N; i++) {
        average_value_arr += arr[i];
    }
    average_value_arr = (1.0 / N) * average_value_arr;

    for (int i = 0; i < N; i++) {
        standard_deviation += pow(arr[i] - average_value_arr, 2);
    }
    standard_deviation = standard_deviation * (1.0 / (N - 1.0));

    double* average_value_arr_i = new double[N];
    double* standard_deviation_i = new double[N];
    double* y = new double[N];
    double average_value_y = 0;
    for (int i = 0; i < N; i++) {
        average_value_arr_i[i] = (1.0 / (N - 1.0)) * (N * average_value_arr - arr[i]);
        standard_deviation_i[i] = (1.0 / (N - 2.0)) * ((N - 1.0) * standard_deviation - (N / (N - 1.0)) * pow(arr[i] - average_value_arr, 2));
        y[i] = pow(standard_deviation_i[i], 1.0 / 3.0);
        average_value_y += y[i];
    }
    average_value_y = (1.0 / N) * average_value_y;
    /*
    r=r1/(r2*r3)^0.5
    */
    double r = 0, r1 = 0, r2 = 0, r3 = 0;
    for (int i = 0; i < N; i++) {
        r1 += (average_value_arr_i[i] - average_value_arr) * (y[i] - average_value_y);
        r2 += pow(average_value_arr_i[i] - average_value_arr, 2);
        r3 += pow(y[i] - average_value_y, 2);
    }
    r = r1 / pow(r2 * r3, 1.0 / 2.0);

    double z = 0, dispersion = 0, coefficient_kurtosis = 0;
    z = (1.0 / 2.0) * log((1 + r) / (1 - r));
    dispersion = (3.0 / N) - (7.324 / pow(N, 2)) + (53.005 / pow(N, 3));
    coefficient_kurtosis = -(11.7 / N) + (55.06 / pow(N, 2));
    /*
    * u(0.95) = 1.96
    * alfa=0.95
    * z_095=z(0.95)
    */
    double z_095 = 0;
    z_095 = sqrt(dispersion) * (1.96 + (1.0 / 24.0) * (pow(1.96, 3) - 3 * 1.96) * coefficient_kurtosis);
    cout << "Критерий Лина-Мудхолкара" << endl;
    if (abs(z) < z_095) {
        cout << "|z|=" << abs(z) << " < z(кр)=" << z_095 << ", нулевая гипотеза нормальности принимается" << endl;
    }
    else {
        cout << "|z|=" << abs(z) << " > z(кр)=" << z_095 << ", нулевая гипотеза нормальности отклоняется" << endl;
    }
    delete[] average_value_arr_i, standard_deviation_i, y;
}

void criterionMartinsIglevich(double* arr, int N) {
    double median_arr = 0.0, median_z = 0.0; 
    double* z_med = new double[N];
    if (N % 2 == 0) {
        median_arr = (arr[N / 2] + arr[(N / 2) - 1]) / 2.0;
    }
    else {
        median_arr = arr[(N / 2) - 1];
    }

    for (int i = 0; i < N; i++) {
        z_med[i] = abs(arr[i] - median_arr);
    }

    BubbleSort(z_med, N);

    if (N % 2 == 0) {
        median_z = (z_med[N / 2] + z_med[(N / 2) - 1]) / 2;
    }
    else {
        median_z = z_med[(N / 2) - 1];
    }

    double* z = new double[N];
    for (int i = 0; i < N; i++) {
        z[i] = (arr[i] - median_arr) / (9 * median_z);
    }
    /*
    * I=(i1*i2^2)/i3
    * I095= I(кр)
    */
    double  i1 = 0, i2 = 0, i3 = 0, I = 0, I095=0;
    for (int i = 0; i < N; i++) {
        i1 += pow(arr[i] - median_arr, 2);
        i2 += (1.0 - pow(z[i], 2)) * (1.0 - 5.0 * pow(z[i], 2));
        i3 += pow(arr[i] - median_arr, 2) * pow((1.0 - pow(z[i],2)),4);   
    }
    I = (1.0 / (N * (N - 1.0))) * ((i1 * pow(i2, 2)) / i3);
    double n = log(N - 1.0);
    if (N > 50) {
        I095 = abs(1.9065 - 2.5465 * n + 0.5652 * pow(n, 2));
    }
    if (N <= 50) {
        I095 = abs(0.7824 - 1.1021 * n + 0.1021 * pow(n, 2));
    }
    cout << "Критерий Мартинеса-Иглевича" << endl;
    if (I < I095) {
        cout << "I=" << I << " < I(кр)=" << I095 << ", нулевая гипотеза нормальности принимается" << endl;
    }
    else {
        cout << "I=" << I << " > I(кр)=" << I095 << ", нулевая гипотеза нормальности отклоняется" << endl;
    }
    delete[] z_med, z;
}

int main(){
    setlocale(0, "");
   /* double time_series[50] = {76.313,77.230,70.510,66.692,65.668,65.312,64.342,64.929,64.601,62.681,64.366,62.201,59.958,58.400,58.109,56.432,57.172,57.831,59.671,
    59.650,57.695,57.731,58.921,58.589,56.788,56.812,57.034,60.462,62.209,62.714,62.883,66.123,67.660,65.887,66.241,67.311,67.347,65.861,
    65.148,64.619,64.816,64.231,63.199,65.533,64.987,64.356,63.865,62.941,61.782,63.884};*/
   double time_series[50] = {83.087,85.910,78.245,75.588,74.268,73.338,71.236,72.790,72.443,69.165,69.640,65.624,63.668,62.176,62.053,60.423,63.097,64.838,68.644,70.396,
   68.804,67.873,69.112,69.361,68.988,70.318,70.355,74.272,73.755,73.224,73.408,76.270,78.963,75.754,75.336,76.653,76.943,74.778,73.755,72.613,
   72.514,72.436,70.990,72.886,71.616,71.065,70.664,69.900,68.725,69.700};
        
    int N = sizeof(time_series) / sizeof(time_series[0]);
    double* time_series_sort = new double[N];
    for (int i = 0; i < N; i++) {
        time_series_sort[i] = time_series[i];
    }
    BubbleSort(time_series_sort, N);
    for (int i = 0; i < N; i++) {
        cout <<setw(6)<<fixed<< setprecision(3)<< time_series[i] << " " << time_series_sort[i] << endl;
    }
    cout << endl;
    criterionShapiroWilk(time_series_sort, N);
    cout << endl;
    criterionLinMudholkar(time_series_sort, N);
    cout << endl;
    criterionMartinsIglevich(time_series_sort, N);

    return 0;
}