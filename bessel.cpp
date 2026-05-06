#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <gsl/gsl_sf_bessel.h>
#include <map>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using cdouble = complex<double>;

// Структура для хранения статистики
struct Statistics {
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    double sum_abs_error = 0.0;
    double sum_rel_error = 0.0;
    int count = 0;
    double x_max_abs = 0.0;
    double x_max_rel = 0.0;
};

// Структура для хранения табличных значений
struct BesselTableEntry {
    double x;
    double J0, J1, J2;
};

// Структура для хранения обоих наборов коэффициентов
struct ChebyshevCoeffs {
    vector<cdouble> a_coeffs;  
    vector<cdouble> b_coeffs;  
};


double map_to_ab(double t, double a, double b) {
    return 0.5 * (b + a) + 0.5 * (b - a) * t;
}

// Загрузка таблиц Абрамовица
map<double, BesselTableEntry> load_abramowitz_table(const string& filename) {
    map<double, BesselTableEntry> table;
    ifstream in(filename);
    
    if (!in.is_open()) {
        cerr << "[!] Не удалось открыть файл " << filename << "\n";
        return table;
    }
    
    string line;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        for (char& c : line) if (c == ',') c = ' ';
        
        istringstream iss(line);
        BesselTableEntry entry;
        
        if (iss >> entry.x >> entry.J0 >> entry.J1 >> entry.J2) {
            table[entry.x] = entry;
        }
    }
    in.close();
    
    cout << "[+] Загружено " << table.size() << " записей из таблиц Абрамовица-Стигана\n";
    return table;
}


// Вычисление функции Ханкеля
cdouble hankel_H1(int nu, double x) {
    if (x <= 0.0) {
        return cdouble(1e30, 1e30);
    }
    double J = gsl_sf_bessel_Jn(nu, x);
    double Y = gsl_sf_bessel_Yn(nu, x);
    return cdouble(J, Y);
}


vector<cdouble> compute_T_coeffs(int nu, double a, double b, int N) {
    vector<cdouble> coeff_T(N, cdouble(0.0, 0.0));
    const double pi = M_PI;

    for (int j = 0; j < N; ++j) {
        double theta = (2.0 * j + 1.0) * pi / (2.0 * N);
        double t_j   = cos(theta);
        double x_j   = map_to_ab(t_j, a, b);
        
        cdouble f_j = hankel_H1(nu, x_j);

        for (int k = 0; k < N; ++k) {
            coeff_T[k] += f_j * cos(k * theta);
        }
    }

    double scale = 2.0 / N;
    for (auto& val : coeff_T) val *= scale;
    coeff_T[0] *= 0.5;  

    return coeff_T;
}


vector<cdouble> compute_U_coeffs(int nu, double a, double b, int N) {
    vector<cdouble> coeff_U(N, cdouble(0.0, 0.0));
    const double pi = M_PI;

    for (int j = 0; j < N; ++j) {
        double theta = (j + 1.0) * pi / (N + 1.0);
        double t_j   = cos(theta);
        double x_j   = map_to_ab(t_j, a, b);
        
        cdouble f_j = hankel_H1(nu, x_j);
        double sin_theta = sin(theta);
        double weight = sin_theta * sin_theta;

        for (int k = 0; k < N; ++k) {
            double U_k = sin((k + 1.0) * theta) / sin_theta;
            coeff_U[k] += f_j * U_k * weight;
        }
    }

    double scale = 2.0 / (N + 1.0);
    for (auto& val : coeff_U) val *= scale;

    return coeff_U;
}


cdouble clenshaw_mixed(const vector<cdouble>& a_coeffs, const vector<cdouble>& b_coeffs, double t) {
    int n_a = a_coeffs.size() - 1;
    cdouble b2_T = 0.0, b1_T = 0.0;
    for (int k = n_a; k > 0; --k) {
        cdouble b0 = a_coeffs[k] + 2.0 * t * b1_T - b2_T;
        b2_T = b1_T;
        b1_T = b0;
    }
    cdouble sum_T = a_coeffs[0] + t * b1_T - b2_T;
    
    int n_b = b_coeffs.size() - 1;
    cdouble b2_U = 0.0, b1_U = 0.0;
    for (int k = n_b; k > 0; --k) {
        cdouble b0 = b_coeffs[k] + 2.0 * t * b1_U - b2_U;
        b2_U = b1_U;
        b1_U = b0;
    }
    cdouble sum_U = b_coeffs[0] + 2.0 * t * b1_U - b2_U;
    
    return 0.5*(sum_T + sum_U);
}


ChebyshevCoeffs compute_mixed_coeffs(int nu, double a, double b, int N) {
    ChebyshevCoeffs coeffs;
    coeffs.a_coeffs = compute_T_coeffs(nu, a, b, N);
    coeffs.b_coeffs = compute_U_coeffs(nu, a, b, N);
    return coeffs;
}


// Сравнение с таблицами и сбор статистики
Statistics compare_with_tables_detailed(
    const map<double, BesselTableEntry>& ref_table,
    int nu, double a, double b, const ChebyshevCoeffs& coeffs,
    bool print_details = false) {
    
    Statistics stats;
    
    if (ref_table.empty()) {
        cout << "[!] Таблицы не загружены.\n";
        return stats;
    }
    
    if (print_details) {
        cout << "\n" << string(90, '=') << "\n";
        cout << "=== ПОДРОБНОЕ СРАВНЕНИЕ ДЛЯ ν = " << nu 
             << " (T_n + U_n) ===\n";
        cout << "Интервал: [" << a << ", " << b << "]\n";
        cout << string(90, '-') << "\n";
        cout << left << setw(8) << "x" 
             << setw(20) << "Таблица" 
             << setw(20) << "Вычислено" 
             << setw(15) << "Абс. погр." 
             << setw(15) << "Отн. погр." << "\n";
        cout << string(90, '-') << "\n";
    }
    
    for (const auto& [x, entry] : ref_table) {
        if (x < a || x > b) continue;
        
        double ref_val = (nu == 0) ? entry.J0 : 
                        (nu == 1) ? entry.J1 : 
                        (nu == 2) ? entry.J2 : 0.0;
        
        if (nu > 2) continue;
        
        double t = 2.0 * (x - a) / (b - a) - 1.0;
        cdouble approx = clenshaw_mixed(coeffs.a_coeffs, coeffs.b_coeffs, t);
        double approx_real = approx.real();
        
        double abs_err = abs(approx_real - ref_val);
        double rel_err = (abs(ref_val) > 1e-10) ? abs_err / abs(ref_val) : abs_err;
        
        stats.max_abs_error = max(stats.max_abs_error, abs_err);
        stats.max_rel_error = max(stats.max_rel_error, rel_err);
        stats.sum_abs_error += abs_err;
        stats.sum_rel_error += rel_err;
        stats.count++;
        
        if (abs_err == stats.max_abs_error) stats.x_max_abs = x;
        if (rel_err == stats.max_rel_error) stats.x_max_rel = x;
        
        if (print_details && (stats.count <= 15 || stats.count % 10 == 0)) {
            cout << fixed << setprecision(2) << setw(8) << x
                 << scientific << setprecision(6)
                 << setw(18) << ref_val
                 << setw(18) << approx_real
                 << setw(13) << abs_err
                 << setw(13) << rel_err << "\n";
        }
    }
    
    return stats;
}


void print_statistics(const Statistics& stats, int nu) {
    cout << "\n" << string(90, '=') << "\n";
    cout << "=== ИТОГОВАЯ СТАТИСТИКА ДЛЯ ν = " << nu << " ===\n";
    cout << string(90, '-') << "\n";
    cout << "Количество точек сравнения:     " << stats.count << "\n";
    cout << "Макс. абсолютная погрешность:   " << scientific << setprecision(6) 
         << stats.max_abs_error << " (при x = " << fixed << setprecision(1) 
         << stats.x_max_abs << ")\n";
    cout << "Макс. относительная погрешность: " << scientific << setprecision(6) 
         << stats.max_rel_error << " (при x = " << fixed << setprecision(1) 
         << stats.x_max_rel << ")\n";
    
    if (stats.count > 0) {
        cout << "Средняя абсолютная погрешность: " << scientific << setprecision(6) 
             << (stats.sum_abs_error / stats.count) << "\n";
        cout << "Средняя относительная погрешность: " << scientific << setprecision(6) 
             << (stats.sum_rel_error / stats.count) << "\n";
    }
    cout << string(90, '=') << "\n";
}


void print_summary_table(const vector<Statistics>& all_stats, const vector<int>& nu_list) {
    cout << "\n\n" << string(100, '=') << "\n";
    cout << "=== СВОДНАЯ ТАБЛИЦА ПОГРЕШНОСТЕЙ (T_n + U_n) ===\n";
    cout << string(100, '=') << "\n";
    cout << left << setw(10) << "ν" 
         << setw(12) << "Точек" 
         << setw(20) << "Макс. |Δ|" 
         << setw(20) << "Макс. Δ/|ref|" 
         << setw(20) << "Сред. |Δ|" 
         << setw(20) << "Сред. Δ/|ref|" << "\n";
    cout << string(100, '-') << "\n";
    
    for (size_t i = 0; i < nu_list.size(); ++i) {
        const auto& stats = all_stats[i];
        cout << scientific << setprecision(6);
        cout << setw(10) << nu_list[i]
             << setw(12) << stats.count
             << setw(18) << stats.max_abs_error
             << setw(18) << stats.max_rel_error
             << setw(18) << (stats.count > 0 ? stats.sum_abs_error/stats.count : 0)
             << setw(18) << (stats.count > 0 ? stats.sum_rel_error/stats.count : 0) << "\n";
    }
    cout << string(100, '=') << "\n";
}


void save_results(const vector<int>& nu_list, 
                  const vector<ChebyshevCoeffs>& all_coeffs,
                  double a, double b, int N) {
    ofstream out("hankel_cheb_mixed_coeffs.csv");
    out << "# Hankel H1 Mixed Chebyshev coefficients (T_n + U_n)\n";
    out << "# Interval: [" << a << ", " << b << "], N = " << N << "\n";
    out << "# Format: nu, k, a_k (for T_k), b_k (for U_k)\n";
    
    for (size_t idx = 0; idx < nu_list.size(); ++idx) {
        int nu = nu_list[idx];
        for (size_t k = 0; k < all_coeffs[idx].a_coeffs.size(); ++k) {
            out << nu << ", " << k << ", " 
                << scientific << setprecision(10)
                << all_coeffs[idx].a_coeffs[k].real() << ", "
                << all_coeffs[idx].a_coeffs[k].imag() << ", "
                << all_coeffs[idx].b_coeffs[k].real() << ", "
                << all_coeffs[idx].b_coeffs[k].imag() << "\n";
        }
    }
    out.close();
    cout << "\n[+] Коэффициенты сохранены в hankel_cheb_mixed_coeffs.csv\n";
}


int main() {
    const int N = 100;
    const double a = 0.5, b = 10.0;
    const int nu_min = 0, nu_max = 2;
    
    cout << "=== ГИБРИДНОЕ РАЗЛОЖЕНИЕ: f ≈ Σa_n*T_n + Σb_n*U_n ===\n";
    cout << "Интервал: [" << a << ", " << b << "], узлов: " << N << "\n";
    cout << "Порядки ν: " << nu_min << " … " << nu_max << "\n";
    cout << "Параллелизм: OpenMP, потоков: " << omp_get_max_threads() << "\n\n";

    string ref_file = "abramowitz_stegan_table.csv";
    auto ref_table = load_abramowitz_table(ref_file);
    
    vector<int> nu_list;
    for (int nu = nu_min; nu <= nu_max; ++nu) nu_list.push_back(nu);
    
    vector<ChebyshevCoeffs> all_coeffs(nu_list.size());
    vector<Statistics> all_stats(nu_list.size());
    
    double start_time = omp_get_wtime();

    #pragma omp parallel for schedule(dynamic)
    for (size_t idx = 0; idx < nu_list.size(); ++idx) {
        int nu = nu_list[idx];
        cout << "[Thread " << omp_get_thread_num() << "] Вычисление a_n и b_n для ν = " << nu << "\n";
        all_coeffs[idx] = compute_mixed_coeffs(nu, a, b, N);
    }

    double elapsed = omp_get_wtime() - start_time;
    cout << "\n[+] Все коэффициенты вычислены за " << fixed << setprecision(3) 
         << elapsed << " с.\n";

    save_results(nu_list, all_coeffs, a, b, N);

    cout << "\n\n";
    for (size_t idx = 0; idx < nu_list.size(); ++idx) {
        all_stats[idx] = compare_with_tables_detailed(
            ref_table, nu_list[idx], a, b, all_coeffs[idx], 
            true
        );
        print_statistics(all_stats[idx], nu_list[idx]);
    }

    print_summary_table(all_stats, nu_list);

    cout << "\n[✓] Тестирование завершено.\n";

    return 0;
}