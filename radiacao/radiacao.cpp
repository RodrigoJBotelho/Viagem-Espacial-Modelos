// radiacao.cpp
// Modelo DIDÁTICO de dose de radiação acumulada em missão tripulada
// em função da duração e da espessura de blindagem.
//
// Modelo (por material):
//   dose_GCR_diaria(x) = D0 * ( f + (1 - f) * exp(-x / lambda_gcr) )   [mSv/dia]
//   dose_SPE(x)        = n_eventos * D_spe * exp(-x / lambda_spe)      [mSv]
//   dose_total(x, t)   = dose_GCR_diaria(x) * t + dose_SPE(x)
//
// x = espessura em densidade superficial (g/cm^2).
// f = fração "residual" da dose que a blindagem não consegue remover
//     (a radiação cósmica galáctica tem componentes muito penetrantes).
//
// ATENÇÃO: os parâmetros abaixo são ILUSTRATIVOS. Substitua por valores
// obtidos de referências abertas (veja README.md) antes de tirar conclusões.
// Este código NÃO substitui simulações de transporte de partículas
// (por exemplo Geant4 ou HZETRN).
//
// Uso: ./radiacao [dias_missao=180] [eventos_solares=1] [espessura_ref=10]
// Saídas (no diretório atual):
//   radiacao_espessura.csv  -> dose total x espessura, para a missão dada
//   radiacao_duracao.csv    -> dose total x duração, para a espessura_ref

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct Material {
    std::string nome;
    double lambda_gcr;  // g/cm^2
    double f_residual;  // adimensional, 0..1
    double lambda_spe;  // g/cm^2
};

// Parâmetros ILUSTRATIVOS (ajustar com dados de literatura).
static const double D0_MSV_DIA = 1.8;   // ordem de grandeza da dose diária em cruzeiro
static const double D_SPE_MSV = 100.0;  // dose por evento solar sem blindagem (muito variável)

static std::vector<Material> materiais() {
    std::vector<Material> m;
    Material al = {"aluminio", 20.0, 0.97, 6.0};
    Material pe = {"polietileno", 15.0, 0.75, 4.0};
    Material ag = {"agua", 15.0, 0.80, 4.5};
    m.push_back(al);
    m.push_back(pe);
    m.push_back(ag);
    return m;
}

static double dose_gcr_diaria(const Material& m, double x) {
    return D0_MSV_DIA * (m.f_residual + (1.0 - m.f_residual) * std::exp(-x / m.lambda_gcr));
}

static double dose_spe(const Material& m, double x, double n_eventos) {
    return n_eventos * D_SPE_MSV * std::exp(-x / m.lambda_spe);
}

static double dose_total(const Material& m, double x, double dias, double n_eventos) {
    return dose_gcr_diaria(m, x) * dias + dose_spe(m, x, n_eventos);
}

int main(int argc, char** argv) {
    double dias = (argc > 1) ? std::atof(argv[1]) : 180.0;
    double n_ev = (argc > 2) ? std::atof(argv[2]) : 1.0;
    double x_ref = (argc > 3) ? std::atof(argv[3]) : 10.0;
    if (dias <= 0 || n_ev < 0 || x_ref < 0) {
        std::cerr << "Uso: ./radiacao [dias>0] [eventos>=0] [espessura_ref>=0]\n";
        return 1;
    }
    std::vector<Material> mats = materiais();

    // 1) dose total x espessura
    {
        std::ofstream f("radiacao_espessura.csv");
        f << "espessura_g_cm2";
        for (size_t j = 0; j < mats.size(); ++j) f << "," << mats[j].nome;
        f << "\n";
        for (int i = 0; i <= 40; ++i) {
            double x = i;
            f << x;
            for (size_t j = 0; j < mats.size(); ++j) f << "," << dose_total(mats[j], x, dias, n_ev);
            f << "\n";
        }
    }

    // 2) dose total x duração (espessura fixa) + sem blindagem
    {
        std::ofstream f("radiacao_duracao.csv");
        f << "dias,sem_blindagem";
        for (size_t j = 0; j < mats.size(); ++j) f << "," << mats[j].nome;
        f << "\n";
        for (int i = 0; i <= 30; ++i) {
            double t = 30.0 * i;
            double ev = (t > 0.0) ? n_ev : 0.0;
            f << t << "," << dose_total(mats[0], 0.0, t, ev);
            for (size_t j = 0; j < mats.size(); ++j) f << "," << dose_total(mats[j], x_ref, t, ev);
            f << "\n";
        }
    }

    std::printf("Missao: %.0f dias, %.1f evento(s) solar(es)\n", dias, n_ev);
    std::printf("%-14s %9s %12s\n", "material", "x[g/cm2]", "dose[mSv]");
    const double xs[] = {0.0, 5.0, 10.0, 20.0};
    for (size_t j = 0; j < mats.size(); ++j)
        for (int i = 0; i < 4; ++i)
            std::printf("%-14s %9.1f %12.1f\n", mats[j].nome.c_str(), xs[i],
                        dose_total(mats[j], xs[i], dias, n_ev));
    std::printf("Arquivos gerados: radiacao_espessura.csv, radiacao_duracao.csv\n");
    return 0;
}
