// NucleonModel.cpp is part of the IP-Glasma solver.

#include "NucleonModel.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "PhysConst.h"

using PhysConst::hbarc;

namespace {
/**
 * Samples from a log-normal distribution with the given mean and
 * variance.
 * \param[in,out] random Random-number source.
 * \param[in] mean Mean of the distribution.
 * \param[in] variance Variance of the distribution.
 * \return The sampled value.
 */
double sampleLogNormalDistribution(
    Random &random, const double mean, const double variance) {
    const double meansq = mean * mean;
    const double mu = log(meansq / sqrt(variance + meansq));
    const double sigma = sqrt(log(variance / meansq + 1.));
    double sampleX = exp(mu + sigma * random.gauss());
    return (sampleX);
}

/**
 * Samples the Qs normalization factors of a nucleon's \p Nq parts: 1
 * each, or log-normal with mean 1 and width \p width if \p smear.
 * \param[in,out] random Random-number source.
 * \param[in] smear Whether the normalization fluctuates.
 * \param[in] width Width of the log-normal fluctuations.
 * \param[in] Nq Number of factors.
 * \return The factors.
 */
std::vector<double> sampleQsNormalization(
    Random &random, bool smear, double width, const int Nq) {
    const double QsSmearWidth = width;
    std::vector<double> gauss_array(Nq, 1.);  // default norm = 1
    if (smear) {
        // introduce a log-normal distribution for Qs normalization
        // dividing by exp(0.5 sigma^2) to ensure the mean is 1
        // the varirance in this case is exp(sigma) - 1 for the log-normal
        // distribution
        for (int iq = 0; iq < Nq; iq++) {
            gauss_array[iq] =
                (exp(random.gauss(0, QsSmearWidth))
                 / exp(QsSmearWidth * QsSmearWidth / 2.));
        }
    }
    return gauss_array;
}
}  // namespace

double GaussianProfile::thickness(double x, double y) const {
    const double xm = xm_;
    const double ym = ym_;
    const double xi = xi_;
    const double phi = phi_;
    double bp2 = (xm - x) * (xm - x) + (ym - y) * (ym - y)
                 + xi * pow((xm - x) * cos(phi) + (ym - y) * sin(phi), 2.);
    bp2 /= hbarc * hbarc;
    return sqrt(1 + xi) * exp(-bp2 / (2. * BG_)) / (2. * M_PI * BG_)
           * normalization_;
}

double HotSpotProfile::thickness(double x, double y) const {
    const double xm = xm_;
    const double ym = ym_;
    double T = 0.;
    for (unsigned int iq = 0; iq < xq_.size(); iq++) {
        double bp2 = (xm + xq_[iq] - x) * (xm + xq_[iq] - x)
                     + (ym + yq_[iq] - y) * (ym + yq_[iq] - y);
        bp2 /= hbarc * hbarc;

        T += exp(-bp2 / (2. * BGq_[iq])) / (2. * M_PI * BGq_[iq])
             / (static_cast<double>(xq_.size())) * normalization_[iq];
    }
    return T;
}

std::unique_ptr<NucleonModel> NucleonModel::create(const Parameters &param) {
    const std::string &model = param.subnucleon.nucleonModel;
    if (model == "gaussian") return std::make_unique<GaussianNucleon>(param);
    if (model == "hotspots") return std::make_unique<HotSpotNucleon>(param);
    throw std::invalid_argument("unknown nucleonModel " + model);
}

GaussianNucleon::GaussianNucleon(const Parameters &param)
    : BG_(param.subnucleon.BG),
      anisotropy_(param.subnucleon.protonAnisotropy),
      smearQs_(param.subnucleon.smearQs),
      smearingWidth_(param.subnucleon.smearingWidth) {}

std::unique_ptr<NucleonProfile> GaussianNucleon::sample(
    Random &random, const ReturnValue &nucleon) const {
    const std::vector<double> normalization =
        sampleQsNormalization(random, smearQs_, smearingWidth_, 1);
    return std::make_unique<GaussianProfile>(
        nucleon.x, nucleon.y, BG_, anisotropy_, nucleon.phi, normalization[0]);
}

HotSpotNucleon::HotSpotNucleon(const Parameters &param)
    : BG_(param.subnucleon.BG),
      BGq_(param.subnucleon.BGq),
      BGqVar_(param.subnucleon.BGqVar),
      dqMin_(param.subnucleon.dqMin),
      omega_(param.subnucleon.omega),
      NqBase_(param.subnucleon.NqBase),
      NqFluc_(param.subnucleon.NqFluc),
      shiftOrigin_(param.subnucleon.shiftConstituentQuarkProtonOrigin),
      smearQs_(param.subnucleon.smearQs),
      smearingWidth_(param.subnucleon.smearingWidth) {}

int HotSpotNucleon::sampleNumberOfPartons(Random &random) const {
    double NqBase = NqBase_;
    int NqBaseInt = static_cast<int>(NqBase);
    double ran = random.genrand64_real2();
    int Nq = NqBaseInt;
    if (ran < NqBase - NqBaseInt) {
        Nq += 1;
    }
    Nq += random.poisson(NqFluc_);
    return (std::max(1, Nq));
}

std::unique_ptr<NucleonProfile> HotSpotNucleon::sample(
    Random &random, const ReturnValue &nucleon) const {
    HotSpotConfiguration hotSpots = sampleHotSpots(random);
    std::vector<double> BGq(hotSpots.x.size(), hotSpots.BGq);
    return std::make_unique<HotSpotProfile>(
        nucleon.x, nucleon.y, std::move(hotSpots.x), std::move(hotSpots.y),
        std::move(BGq), std::move(hotSpots.normalization));
}

HotSpotConfiguration HotSpotNucleon::sampleHotSpots(
    Random &random, int number) const {
    std::vector<double> x_array, y_array, z_array;
    const double sqrtBG = sqrt(BG_) * hbarc;  // fm
    const double BGqMean = BGq_;
    const double BGqVar = BGqVar_;
    const double BGq =
        (0.09 + sampleLogNormalDistribution(random, BGqMean - 0.09, BGqVar));
    const int Nq = (number > 0) ? number : sampleNumberOfPartons(random);
    const double dq_min = dqMin_;  // fm
    const double dq_min_sq = dq_min * dq_min;
    const double omega = omega_;

    std::vector<double> r_array(Nq, 0.);
    for (int iq = 0; iq < Nq; iq++) {
        if (std::abs(omega - 1) < 1e-8) {
            double xq = sqrtBG * random.gauss();
            double yq = sqrtBG * random.gauss();
            double zq = sqrtBG * random.gauss();
            r_array[iq] = sqrt(xq * xq + yq * yq + zq * zq);
        } else {
            double bperp = sqrtBG * sqrt(omega * random.sampleGammaInc());
            r_array[iq] = bperp;  // bperp in 2D (asuume z = 0)
        }
    }
    std::sort(r_array.begin(), r_array.end());

    x_array.resize(Nq, 0.);
    y_array.resize(Nq, 0.);
    z_array.resize(Nq, 0.);
    for (unsigned int i = 0; i < r_array.size(); i++) {
        double r_i = r_array[i];
        int reject_flag = 0;
        int iter = 0;
        double x_i, y_i, z_i;
        do {
            iter++;
            reject_flag = 0;
            double phi = 2. * M_PI * random.genrand64_real2();
            double theta = acos(1. - 2. * random.genrand64_real2());
            if (std::abs(omega - 1) < 1e-8) {
                x_i = r_i * sin(theta) * cos(phi);
                y_i = r_i * sin(theta) * sin(phi);
                z_i = r_i * cos(theta);
            } else {
                x_i = r_i * cos(phi);
                y_i = r_i * sin(phi);
                z_i = 0.;  // assume z=0
            }
            for (int j = i - 1; j >= 0; j--) {
                if ((r_i - r_array[j]) * (r_i - r_array[j]) > dq_min_sq) break;
                double dsq =
                    ((x_i - x_array[j]) * (x_i - x_array[j])
                     + (y_i - y_array[j]) * (y_i - y_array[j])
                     + (z_i - z_array[j]) * (z_i - z_array[j]));
                if (dsq < dq_min_sq) {
                    reject_flag = 1;
                    break;
                }
            }
        } while (reject_flag == 1 && iter < 100);
        x_array[i] = x_i;
        y_array[i] = y_i;
        z_array[i] = z_i;
    }
    double avgxq = 0.;
    double avgyq = 0.;
    double avgzq = 0.;
    if (shiftOrigin_) {
        for (int iq = 0; iq < Nq; iq++) {
            avgxq += x_array[iq];
            avgyq += y_array[iq];
            avgzq += z_array[iq];
        }
        avgxq /= static_cast<double>(Nq);
        avgyq /= static_cast<double>(Nq);
        avgzq /= static_cast<double>(Nq);
        for (int iq = 0; iq < Nq; iq++) {
            x_array[iq] -= avgxq;
            y_array[iq] -= avgyq;
            z_array[iq] -= avgzq;
        }
    }

    const int Npartons = std::max(1, static_cast<int>(x_array.size()));
    HotSpotConfiguration hotSpots;
    hotSpots.normalization =
        sampleQsNormalization(random, smearQs_, smearingWidth_, Npartons);
    hotSpots.x = std::move(x_array);
    hotSpots.y = std::move(y_array);
    hotSpots.z = std::move(z_array);
    hotSpots.BGq = BGq;
    return hotSpots;
}
