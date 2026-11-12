// calor1d.cpp
// Condução de calor 1D transiente em uma placa (modelo simplificado de escudo
// térmico de reentrada):   dT/dt = alpha * d2T/dx2
//
// Contorno:  x = 0 (face externa)  -> fluxo de calor imposto  -k dT/dx = -q
//            x = L (face interna)  -> isolada (dT/dx = 0)
//
// Métodos:   (1) FTCS explícito  (estável se r = alpha*dt/dx^2 <= 0,5)
//            (2) Crank-Nicolson implícito (incondicionalmente estável)
//
// Validação: solução analítica do meio semi-infinito com fluxo constante
//   T(x,t) = T0 + (2q/k) sqrt(alpha t/pi) exp(-x^2/(4 alpha t))
//                - (q x / k) erfc( x / (2 sqrt(alpha t)) )
// Vale enquanto a frente térmica não alcança x = L (aqui L >> 2 sqrt(alpha tf)).
//
// Parâmetros de material ILUSTRATIVOS (k, alpha, q). Substitua por dados reais.
//
// Uso: ./calor1d [N=100] [dt=0.4*dx^2/alpha]
// Saídas (no diretório atual): reentrada_perfil.csv, reentrada_historico.csv

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

typedef std::vector<double> Vec;
static const double PI = 3.14159265358979323846;

struct Params {
    double L, k, alpha, q, T0, tf;
    Params()
        : L(0.10),       // espessura da placa [m]
          k(1.0),        // condutividade [W/(m K)]
          alpha(1e-6),   // difusividade térmica [m^2/s]
          q(5.0e4),      // fluxo de calor na face externa [W/m^2]
          T0(300.0),     // temperatura inicial [K]
          tf(60.0) {}    // tempo final [s]
};

static double analitica(const Params& p, double x, double t) {
    if (t <= 0.0) return p.T0;
    double s = std::sqrt(p.alpha * t);
    return p.T0 + (2.0 * p.q / p.k) * (s / std::sqrt(PI)) * std::exp(-x * x / (4.0 * p.alpha * t))
           - (p.q * x / p.k) * std::erfc(x / (2.0 * s));
}

// Um passo explícito (nó fantasma em x=0 para impor o fluxo; isolado em x=L)
static void passo_explicito(const Params& p, const Vec& T, Vec& Tn, double r, double dx) {
    int N = (int)T.size() - 1;
    Tn[0] = T[0] + r * (2.0 * T[1] - 2.0 * T[0] + 2.0 * dx * p.q / p.k);
    for (int i = 1; i < N; ++i) Tn[i] = T[i] + r * (T[i + 1] - 2.0 * T[i] + T[i - 1]);
    Tn[N] = T[N] + r * (2.0 * T[N - 1] - 2.0 * T[N]);
}

// Algoritmo de Thomas para sistema tridiagonal: a (inferior), b (diagonal), c (superior)
static void thomas(Vec a, Vec b, Vec c, Vec d, Vec& x) {
    int n = (int)b.size();
    for (int i = 1; i < n; ++i) {
        double w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    x.assign(n, 0.0);
    x[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i) x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
}

// Um passo theta-method (theta = 0.5 -> Crank-Nicolson)
static void passo_cn(const Params& p, const Vec& T, Vec& Tn, double r, double dx, double theta) {
    int n = (int)T.size();
    int N = n - 1;
    Vec a(n, -theta * r), b(n, 1.0 + 2.0 * theta * r), c(n, -theta * r), d(n);
    c[0] = -2.0 * theta * r;
    a[N] = -2.0 * theta * r;
    // lado direito: (I + (1-theta) r A) T + r * s, com s0 = 2 dx q / k
    d[0] = T[0] + (1.0 - theta) * r * (-2.0 * T[0] + 2.0 * T[1]) + r * 2.0 * dx * p.q / p.k;
    for (int i = 1; i < N; ++i)
        d[i] = T[i] + (1.0 - theta) * r * (T[i + 1] - 2.0 * T[i] + T[i - 1]);
    d[N] = T[N] + (1.0 - theta) * r * (2.0 * T[N - 1] - 2.0 * T[N]);
    thomas(a, b, c, d, Tn);
}

int main(int argc, char** argv) {
    Params p;
    int N = (argc > 1) ? std::atoi(argv[1]) : 100;
    if (N < 10) { std::cerr << "N deve ser >= 10\n"; return 1; }
    double dx = p.L / N;
    double dt = (argc > 2) ? std::atof(argv[2]) : 0.4 * dx * dx / p.alpha;
    if (dt <= 0.0) { std::cerr << "dt deve ser > 0\n"; return 1; }
    int nsteps = (int)std::ceil(p.tf / dt - 1e-12);
    dt = p.tf / nsteps;
    double r = p.alpha * dt / (dx * dx);

    std::printf("N=%d  dx=%.4g m  dt=%.4g s  passos=%d  r=%.3f\n", N, dx, dt, nsteps, r);
    if (r > 0.5)
        std::printf("AVISO: r > 0,5 -> o metodo explicito e INSTAVEL (o CN continua estavel).\n");

    Vec Te(N + 1, p.T0), Tc(N + 1, p.T0), tmp(N + 1);
    std::ofstream h("reentrada_historico.csv");
    h << "t_s,T_superficie_explicito,T_superficie_crank_nicolson,T_superficie_analitica\n";
    int salto = std::max(1, nsteps / 500);
    double errH_e = 0, errH_c = 0;
    h << 0 << "," << Te[0] << "," << Tc[0] << "," << p.T0 << "\n";
    for (int n = 1; n <= nsteps; ++n) {
        passo_explicito(p, Te, tmp, r, dx);
        Te = tmp;
        passo_cn(p, Tc, tmp, r, dx, 0.5);
        Tc = tmp;
        double t = n * dt, Ta = analitica(p, 0.0, t);
        errH_e = std::max(errH_e, std::fabs(Te[0] - Ta));
        errH_c = std::max(errH_c, std::fabs(Tc[0] - Ta));
        if (n % salto == 0 || n == nsteps)
            h << t << "," << Te[0] << "," << Tc[0] << "," << Ta << "\n";
    }

    std::ofstream f("reentrada_perfil.csv");
    f << "x_m,T_explicito,T_crank_nicolson,T_analitica\n";
    double errP_e = 0, errP_c = 0;
    for (int i = 0; i <= N; ++i) {
        double x = i * dx, Ta = analitica(p, x, p.tf);
        errP_e = std::max(errP_e, std::fabs(Te[i] - Ta));
        errP_c = std::max(errP_c, std::fabs(Tc[i] - Ta));
        f << x << "," << Te[i] << "," << Tc[i] << "," << Ta << "\n";
    }

    std::printf("T na superficie em t=%.0f s: explicito=%.2f  CN=%.2f  analitica=%.2f K\n",
                p.tf, Te[0], Tc[0], analitica(p, 0.0, p.tf));
    std::printf("Erro maximo no perfil final : explicito=%.3g K   CN=%.3g K\n", errP_e, errP_c);
    std::printf("Erro maximo na superficie   : explicito=%.3g K   CN=%.3g K\n", errH_e, errH_c);
    std::printf("Arquivos gerados: reentrada_perfil.csv, reentrada_historico.csv\n");
    return 0;
}
