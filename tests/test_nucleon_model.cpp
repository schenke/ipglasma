#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include "Glauber.h"
#include "Init.h"
#include "NucleonModel.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "Random.h"
#include "doctest.h"

using PhysConst::hbarc;

namespace {
// Integral of a profile over the transverse plane, in units of the
// thickness normalization (d^2b in GeV^-2).
double integrate(const NucleonProfile &profile, double xm, double ym) {
    const double range = 4.;   // fm
    const double step = 0.02;  // fm
    double sum = 0.;
    for (double x = xm - range; x < xm + range; x += step) {
        for (double y = ym - range; y < ym + range; y += step) {
            sum += profile.thickness(x, y);
        }
    }
    return sum * (step / hbarc) * (step / hbarc);
}

ReturnValue nucleonAt(double x, double y, double phi = 0.) {
    ReturnValue nucleon {};
    nucleon.x = x;
    nucleon.y = y;
    nucleon.phi = phi;
    return nucleon;
}

void makeHotSpotParam(Parameters &param) {
    param.subnucleon.nucleonModel = "hotspots";
    param.subnucleon.BG = 4.;
    param.subnucleon.BGq = 0.3;
    param.subnucleon.BGqVar = 0.;
    param.subnucleon.dqMin = 0.;
    param.subnucleon.omega = 1.;
    param.subnucleon.NqBase = 3.;
    param.subnucleon.NqFluc = 0.;
    param.subnucleon.shiftConstituentQuarkProtonOrigin = true;
    param.subnucleon.smearQs = false;
    param.subnucleon.smearingWidth = 0.;
}
}  // namespace

TEST_CASE("GaussianProfile: 1/(2 pi BG) at the center, integrates to 1") {
    const GaussianProfile profile(0.5, -0.3, 1.0, 0., 0., 1.);
    CHECK(profile.thickness(0.5, -0.3) == doctest::Approx(1.0 / (2.0 * M_PI)));
    CHECK(integrate(profile, 0.5, -0.3) == doctest::Approx(1.).epsilon(1e-4));

    // an elongated nucleon keeps its normalization, a Qs factor scales it
    const GaussianProfile elongated(0., 0., 4.0, 0.8, 0.7, 1.3);
    CHECK(integrate(elongated, 0., 0.) == doctest::Approx(1.3).epsilon(1e-4));
}

TEST_CASE("HotSpotProfile: one hot spot at the center is a Gaussian") {
    const HotSpotProfile profile(0.2, 0.1, {0.}, {0.}, {1.0}, {1.});
    CHECK(profile.numberOfHotSpots() == 1);
    CHECK(profile.thickness(0.2, 0.1) == doctest::Approx(1.0 / (2.0 * M_PI)));
}

TEST_CASE(
    "HotSpotNucleon: samples Nq hot spots whose profile integrates to 1") {
    Parameters param;
    makeHotSpotParam(param);
    const HotSpotNucleon model(param);
    Random random;
    random.init_genrand64(7ULL);

    for (int i = 0; i < 5; ++i) {
        const std::unique_ptr<NucleonProfile> profile =
            model.sample(random, nucleonAt(1.0, -0.5));
        const auto *hotSpots =
            dynamic_cast<const HotSpotProfile *>(profile.get());
        REQUIRE(hotSpots != nullptr);
        CHECK(hotSpots->numberOfHotSpots() == 3);
        CHECK(
            integrate(*profile, 1.0, -0.5)
            == doctest::Approx(1.).epsilon(1e-3));
    }
}

TEST_CASE("HotSpotNucleon: a fractional NqBase gives the right mean") {
    Parameters param;
    makeHotSpotParam(param);
    param.subnucleon.NqBase = 2.5;
    const HotSpotNucleon model(param);
    Random random;
    random.init_genrand64(11ULL);

    const int samples = 20000;
    double sum = 0.;
    for (int i = 0; i < samples; ++i) {
        const int n = model.sampleNumberOfPartons(random);
        CHECK((n == 2 || n == 3));
        sum += n;
    }
    CHECK(sum / samples == doctest::Approx(2.5).epsilon(0.02));
}

TEST_CASE("NucleonModel::create returns the selected model") {
    Parameters param;
    makeHotSpotParam(param);
    CHECK(
        dynamic_cast<HotSpotNucleon *>(NucleonModel::create(param).get())
        != nullptr);
    param.subnucleon.nucleonModel = "gaussian";
    CHECK(
        dynamic_cast<GaussianNucleon *>(NucleonModel::create(param).get())
        != nullptr);
    param.subnucleon.nucleonModel = "unknown";
    CHECK_THROWS_AS(NucleonModel::create(param), std::invalid_argument);
}

TEST_CASE(
    "Init::computeNucleonThicknessAtCell sums the nucleons' thickness and "
    "divides by the number of averaged nuclei") {
    std::vector<std::unique_ptr<NucleonProfile>> profiles;
    profiles.push_back(
        std::make_unique<GaussianProfile>(0., 0., 1.0, 0., 0., 1.));
    profiles.push_back(
        std::make_unique<GaussianProfile>(0., 0., 1.0, 0., 0., 2.));
    int nn[2] = {4, 4};
    Init init(nn);
    CHECK(
        init.computeNucleonThicknessAtCell(profiles, 0., 0., 2.)
        == doctest::Approx(3. / (2.0 * M_PI) / 2.));
}
