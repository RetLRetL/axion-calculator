#include <iostream>
#include <cstring> // for sprintf
#include <cmath>

#include "TFile.h"
#include "TH1D.h"

#include "settings.h"
#include "database.h"

using namespace std;

int main()
{
    load_XCOM();
    init_mixture();

    const double gagg = 1e-2;
    const double Z = 92.0;
    const double me = 0.5;
    const double ma = 1.0;
    const double mU = 221700.0;

    const double alpha = 111.0 * pow(Z, -1.0/3.0) * me;
    const double a = alpha;

    // ===== FILE PATHS =====
    char input_file[256];
    char output_file[256];

    sprintf(input_file, "%s/Gamma_sp.root", DATA_DIR);
    sprintf(output_file, "%s/alp_result.root", DATA_DIR);

    TFile *f1 = new TFile(input_file, "READ");
    TFile *fout = new TFile(output_file, "RECREATE");

    TH1D *h1 = (TH1D*)f1->Get("Sp_gamma_start");
    TH1D *h2 = (TH1D*)f1->Get("Sp_gamma_end");

    TH1D *h1_out = (TH1D*)h1->Clone("alp_start");
    TH1D *h2_out = (TH1D*)h2->Clone("alp_end");

    h1_out->Reset();
    h2_out->Reset();

    h1_out->GetXaxis()->SetTitle("E, MeV");
    h1_out->GetYaxis()->SetTitle("dN/dE, MeV^{-1} c^{-1} ton^{-1}");

    h2_out->GetXaxis()->SetTitle("E, MeV");
    h2_out->GetYaxis()->SetTitle("dN/dE, MeV^{-1} c^{-1} ton^{-1}");

    int Nbins = h1->GetNbinsX();

    // ===== CROSS SECTIONS =====
    auto sigma_PC = [&](double Eg){
        return gagg*gagg * alpha * Z*Z / 4.0 *
               (log(4.0 * Eg * Eg / (me*me)) - 1.0);
    };

    auto dsigma_dEa = [&](double Eg, double Ea){
        double arg = 1.0 + (2.0 * Eg * (Eg - Ea)) / (a*a);
        double denom = pow(arg, 4);
        double result = gagg*gagg * alpha * Z*Z / (8.0 * pow(Eg,4)) / denom;

        if (USE_CUTS)
        {
            if (isnan(result) || isinf(result)) result = 0;
            if (result > 1e20) result = 1e20;
        }

        return result;
    };

    // ===== MAIN LOOP =====
    for (int i = 1; i <= Nbins; i++)
    {
        double Ea = h1->GetBinCenter(i);

        // protect against sqrt(Ea^2 - ma^2)
        if (Ea <= ma) continue;

        double sqrt_term = sqrt(Ea*Ea - ma*ma);

        double denom_plus  = (mU - Ea + sqrt_term);
        double denom_minus = (mU - Ea - sqrt_term);

        // avoid division by zero
        if (fabs(denom_plus) < 1e-12 || fabs(denom_minus) < 1e-12)
            continue;

        double Egamma_min_dyn = (-0.5 * ma*ma + Ea * mU) / denom_plus;
        double Egamma_max_dyn = (-0.5 * ma*ma + Ea * mU) / denom_minus;

        if (Egamma_min_dyn > Egamma_max_dyn)
            continue;

        // additional safety checks
        if (Egamma_max_dyn <= 0) continue;
        if (Egamma_min_dyn < 0) Egamma_min_dyn = 0;

        double integral1 = 0.0;
        double integral2 = 0.0;

        int used_bins = 0;

        for (int j = 1; j <= Nbins; j++)
        {
            double Eg = h1->GetBinCenter(j);

            if (Eg < Egamma_min_dyn || Eg > Egamma_max_dyn) continue;

            double dEg = h1->GetBinWidth(j);

            double flux1 = h1->GetBinContent(j);
            double flux2 = h2->GetBinContent(j);

            if (flux1 < MIN_VAL && flux2 < MIN_VAL) continue;

            double sigma_tot = sigma_SM_mix(Eg) + sigma_PC(Eg);
            if (sigma_tot <= 0) continue;

            double kernel = dsigma_dEa(Eg, Ea) / sigma_tot;

            integral1 += kernel * flux1 * dEg;
            integral2 += kernel * flux2 * dEg;

            used_bins++;
        }

        // fallback: if no bins fall into the integration range
        if (used_bins == 0)
        {
            double Eg_mid = 0.5 * (Egamma_min_dyn + Egamma_max_dyn);
            double dEg_eff = Egamma_max_dyn - Egamma_min_dyn;

            int bin = h1->FindBin(Eg_mid);

            if (bin >= 1 && bin <= Nbins)
            {
                double flux1 = h1->GetBinContent(bin);
                double flux2 = h2->GetBinContent(bin);

                if (flux1 >= MIN_VAL || flux2 >= MIN_VAL)
                {
                    double sigma_tot = sigma_SM_mix(Eg_mid) + sigma_PC(Eg_mid);

                    if (sigma_tot > 0)
                    {
                        double kernel = dsigma_dEa(Eg_mid, Ea) / sigma_tot;

                        integral1 += kernel * flux1 * dEg_eff;
                        integral2 += kernel * flux2 * dEg_eff;

                        cout << "[FALLBACK USED] Ea=" << Ea
                             << " Eg_mid=" << Eg_mid
                             << " width=" << dEg_eff << endl;
                    }
                }
            }
        }

        h1_out->SetBinContent(i, integral1);
        h2_out->SetBinContent(i, integral2);
    }

    fout->Write();
    fout->Close();
    f1->Close();

    cout << "Done!" << endl;
}