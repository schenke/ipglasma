#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

#include "Glauber.h"
#include "NucleusSampler.h"
#include "Parameters.h"
#include "Random.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
Nucleus makeNucleus(int A, int Z) {
    Nucleus data {};
    data.A = A;
    data.Z = Z;
    data.R_WS = 1.1 * std::cbrt(static_cast<double>(A));
    data.a_WS = 0.5;
    data.d_min = 0.9;
    return data;
}

int countProtons(const std::vector<ReturnValue> &nucleus) {
    return static_cast<int>(std::count_if(
        nucleus.begin(), nucleus.end(),
        [](const ReturnValue &n) { return n.proton; }));
}

void checkCentered(const std::vector<ReturnValue> &nucleus) {
    double x = 0., y = 0., z = 0.;
    for (const auto &n : nucleus) {
        x += n.x;
        y += n.y;
        z += n.z;
    }
    CHECK(std::abs(x) < 1e-10);
    CHECK(std::abs(y) < 1e-10);
    CHECK(std::abs(z) < 1e-10);
}

double minimumDistance(const std::vector<ReturnValue> &nucleus) {
    double dmin = 1e30;
    for (size_t i = 0; i < nucleus.size(); i++) {
        for (size_t j = 0; j < i; j++) {
            const double dx = nucleus[i].x - nucleus[j].x;
            const double dy = nucleus[i].y - nucleus[j].y;
            const double dz = nucleus[i].z - nucleus[j].z;
            dmin = std::min(dmin, std::sqrt(dx * dx + dy * dy + dz * dz));
        }
    }
    return dmin;
}

void checkSamePositions(
    const std::vector<ReturnValue> &a, const std::vector<ReturnValue> &b) {
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); i++) {
        CHECK(a[i].x == b[i].x);
        CHECK(a[i].y == b[i].y);
        CHECK(a[i].z == b[i].z);
        CHECK(a[i].proton == b[i].proton);
    }
}
}  // namespace

TEST_CASE("NucleusSampler::sample places a proton at the origin") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.collision.nucleiToAverage = 1;
    Glauber glauber;
    glauber.initGlauber(
        4.2, "p", "p", /*inb=*/0.0, /*setWSDeformParams=*/false, 0., 0., 0., 0.,
        0., 0., /*forceDminFlag=*/false, 0., 0., 0., /*imax=*/1000);
    Random random;
    random.init_genrand64(1);

    NucleusSampler sampler;
    const Nuclei nuclei = sampler.sample(&param, &random, &glauber);
    const std::vector<ReturnValue> &nucleusA = nuclei.projectile;
    const std::vector<ReturnValue> &nucleusB = nuclei.target;
    for (const auto *nucleus : {&nucleusA, &nucleusB}) {
        REQUIRE(nucleus->size() == 1);
        CHECK(nucleus->front().x == 0.);
        CHECK(nucleus->front().y == 0.);
        CHECK(nucleus->front().z == 0.);
        CHECK(nucleus->front().proton);
        CHECK(nucleus->front().collided == 0);
    }
}

TEST_CASE(
    "NucleusSampler::sample overlays nucleiToAverage independent nuclei") {
    Parameters param;
    makeInitTestParam(param, 4);
    // a transverse polarization turns every nucleus the same way without
    // random numbers, so the first projectile does not depend on how many
    // nuclei follow it
    param.nucleus.polarizationProjectile = 2;
    param.nucleus.polarizationTarget = 2;
    Glauber glauber;
    glauber.initGlauber(
        42., "O", "Ca", /*inb=*/0.0, /*setWSDeformParams=*/false, 0., 0., 0.,
        0., 0., 0., /*forceDminFlag=*/false, 0., 0., 0., /*imax=*/1000);
    auto sampleWith = [&](int nucleiToAverage) {
        param.collision.nucleiToAverage = nucleiToAverage;
        Random random;
        random.init_genrand64(3);
        NucleusSampler sampler;
        return sampler.sample(&param, &random, &glauber);
    };
    const Nuclei single = sampleWith(1);
    const Nuclei averaged = sampleWith(3);

    struct Kind {
        const std::vector<ReturnValue> &single, &averaged;
        int A, Z;
    };
    for (const Kind &kind :
         {Kind {single.projectile, averaged.projectile, 40, 20},
          Kind {single.target, averaged.target, 16, 8}}) {
        CAPTURE(kind.A);
        REQUIRE(kind.single.size() == static_cast<size_t>(kind.A));
        REQUIRE(kind.averaged.size() == static_cast<size_t>(3 * kind.A));
        std::vector<std::vector<ReturnValue>> nuclei;
        for (int i = 0; i < 3; i++) {
            nuclei.emplace_back(
                kind.averaged.begin() + i * kind.A,
                kind.averaged.begin() + (i + 1) * kind.A);
            CHECK(countProtons(nuclei.back()) == kind.Z);
            checkCentered(nuclei.back());
        }
        // independently sampled
        CHECK(nuclei[0].front().x != nuclei[1].front().x);
        CHECK(nuclei[1].front().x != nuclei[2].front().x);
    }
    checkSamePositions(
        std::vector<ReturnValue>(
            averaged.projectile.begin(), averaged.projectile.begin() + 40),
        single.projectile);
}

TEST_CASE(
    "NucleusSampler::generate gives A centered nucleons with Z protons for "
    "every Woods-Saxon variant") {
    struct Variant {
        const char *name;
        double beta2, gamma;
        bool forceDmin;
    };
    for (const Variant &v :
         {Variant {"spherical", 0., 0., false},
          Variant {"axially deformed", 0.3, 0., false},
          Variant {"triaxial", 0.3, 0.5, false},
          Variant {"forced d_min", 0.3, 0.5, true}}) {
        CAPTURE(v.name);
        Nucleus data = makeNucleus(40, 18);
        data.beta2 = v.beta2;
        data.gamma = v.gamma;
        data.forceDminFlag = v.forceDmin;
        data.dR_np = 0.1;
        data.da_np = 0.05;
        Random random;
        random.init_genrand64(7);
        NucleusSampler sampler;
        const std::vector<ReturnValue> nucleus =
            sampler.generate(&random, data);

        REQUIRE(nucleus.size() == 40);
        CHECK(countProtons(nucleus) == 18);
        checkCentered(nucleus);
        for (const auto &n : nucleus) CHECK(n.collided == 0);
        if (v.forceDmin) CHECK(minimumDistance(nucleus) >= data.d_min);
    }
}

TEST_CASE(
    "NucleusSampler::generate dispatches on the deformation and the "
    "forced-d_min flag") {
    using Generator =
        std::vector<ReturnValue> (NucleusSampler::*)(Random *, const Nucleus &);
    struct Case {
        double beta2, beta3, beta4, gamma;
        bool forceDmin;
        Generator expected;
    };
    for (const Case &c :
         {Case {0., 0., 0., 0., false, &NucleusSampler::generateWoodsSaxon},
          // d_min is not forced without deformation
          Case {0., 0., 0., 0., true, &NucleusSampler::generateWoodsSaxon},
          Case {
              0., 0.1, 0., 0., false,
              &NucleusSampler::generateDeformedWoodsSaxon},
          Case {
              0., 0., 0.1, 0., false,
              &NucleusSampler::generateDeformedWoodsSaxon},
          Case {
              0.2, 0., 0., 0.4, false,
              &NucleusSampler::generateTriaxialWoodsSaxon},
          Case {
              0.2, 0., 0., 0., true,
              &NucleusSampler::generateDeformedWoodsSaxonForceDmin}}) {
        CAPTURE(c.beta2);
        CAPTURE(c.beta3);
        CAPTURE(c.beta4);
        CAPTURE(c.gamma);
        CAPTURE(c.forceDmin);
        Nucleus data = makeNucleus(16, 8);
        data.beta2 = c.beta2;
        data.beta3 = c.beta3;
        data.beta4 = c.beta4;
        data.gamma = c.gamma;
        data.forceDminFlag = c.forceDmin;
        NucleusSampler sampler;
        Random random1, random2;
        random1.init_genrand64(11);
        random2.init_genrand64(11);
        const std::vector<ReturnValue> dispatched =
            sampler.generate(&random1, data);
        const std::vector<ReturnValue> direct =
            (sampler.*c.expected)(&random2, data);
        checkSamePositions(dispatched, direct);
    }
}

TEST_CASE(
    "NucleusSampler::sampleFromConfigurations recenters one tabulated "
    "configuration and falls back to Woods-Saxon without one") {
    // two configurations of 4 nucleons, the second shifted by (10, 20, 30)
    std::vector<std::vector<float>> configs(2);
    for (int c = 0; c < 2; c++) {
        for (int i = 0; i < 4; i++) {
            configs[c].push_back(static_cast<float>(i + 10 * c));
            configs[c].push_back(static_cast<float>(2 * i + 20 * c));
            configs[c].push_back(static_cast<float>(-i + 30 * c));
        }
    }
    const Nucleus data = makeNucleus(4, 2);
    NucleusSampler sampler;
    Random random;
    random.init_genrand64(3);
    const std::vector<ReturnValue> nucleus =
        sampler.sampleFromConfigurations(&random, data, configs);

    REQUIRE(nucleus.size() == 4);
    CHECK(countProtons(nucleus) == 2);
    checkCentered(nucleus);
    // both configurations are (i, 2i, -i) around their center 1.5 (1, 2, -1)
    std::vector<double> xs;
    for (const auto &n : nucleus) {
        CHECK(n.y == doctest::Approx(2. * n.x));
        CHECK(n.z == doctest::Approx(-n.x));
        xs.push_back(n.x);
    }
    std::sort(xs.begin(), xs.end());
    for (int i = 0; i < 4; i++) CHECK(xs[i] == doctest::Approx(i - 1.5));

    const std::vector<ReturnValue> generated =
        sampler.sampleFromConfigurations(&random, data, {});
    CHECK(generated.size() == 4);
    CHECK(countProtons(generated) == 2);
}

TEST_CASE(
    "NucleusSampler::assignProtons labels exactly Z nucleons and keeps the "
    "positions") {
    std::vector<ReturnValue> nucleus(10);
    for (int i = 0; i < 10; i++) {
        nucleus[i] = ReturnValue {};
        nucleus[i].x = i;
    }
    Random random;
    random.init_genrand64(5);
    NucleusSampler::assignProtons(&random, nucleus, 4);
    CHECK(countProtons(nucleus) == 4);
    std::vector<double> xs;
    for (const auto &n : nucleus) xs.push_back(n.x);
    std::sort(xs.begin(), xs.end());
    for (int i = 0; i < 10; i++) CHECK(xs[i] == i);
}

TEST_CASE(
    "NucleusSampler rotations preserve distances; transverse polarization "
    "turns the z axis into +y") {
    std::vector<ReturnValue> axes(3, ReturnValue {});
    axes[0].x = 1.;
    axes[1].y = 1.;
    axes[2].z = 1.;
    std::vector<ReturnValue> rotated = axes;
    NucleusSampler::applyPolarizationRotation(nullptr, 2, rotated);
    CHECK(rotated[2].x == doctest::Approx(0.));
    CHECK(rotated[2].y == doctest::Approx(1.));
    CHECK(rotated[2].z == doctest::Approx(0.));

    Random random;
    random.init_genrand64(9);
    for (int flag : {0, 1}) {
        CAPTURE(flag);
        rotated = axes;
        NucleusSampler::applyPolarizationRotation(&random, flag, rotated);
        for (const auto &n : rotated) {
            CHECK(n.x * n.x + n.y * n.y + n.z * n.z == doctest::Approx(1.));
        }
        CHECK(minimumDistance(rotated) == doctest::Approx(std::sqrt(2.)));
        // longitudinal polarization only rotates around the beam axis
        if (flag == 1) CHECK(rotated[2].z == doctest::Approx(1.));
    }
}

TEST_CASE(
    "NucleusSampler::readConfigurationFile keeps only (x, y, z) from the "
    "4-entry "
    "Au197 configuration format") {
    Parameters param;
    makeInitTestParam(param, 4);
    const std::string dir = "test_nucleus_configs_tmp";
    const std::string file = dir + "/Au197.bin.in";
    REQUIRE(std::system(("mkdir -p " + dir).c_str()) == 0);

    // Two configurations, each 197 nucleons of (x, y, z, flag), with
    // values encoding (config, nucleon, component) so misalignment shows.
    const int A = 197;
    const int nConfigs = 2;
    {
        std::ofstream out(file, std::ios::binary);
        for (int c = 0; c < nConfigs; c++) {
            for (int i = 0; i < A; i++) {
                for (int j = 0; j < 4; j++) {
                    float v = (j == 3)
                                  ? -1.f
                                  : static_cast<float>(1000 * c + 3 * i + j);
                    out.write(reinterpret_cast<const char *>(&v), sizeof(v));
                }
            }
        }
    }
    param.nucleus.nuclearConfigurationsPath = dir;

    NucleusSampler sampler;
    const std::vector<std::vector<float>> configs =
        sampler.readConfigurationFile(A, 0, 0, 0., nullptr, &param);

    REQUIRE(configs.size() == static_cast<size_t>(nConfigs));
    for (int c = 0; c < nConfigs; c++) {
        REQUIRE(configs[c].size() == static_cast<size_t>(3 * A));
        for (int i = 0; i < A; i++) {
            for (int j = 0; j < 3; j++) {
                CHECK(configs[c][3 * i + j] == 1000 * c + 3 * i + j);
            }
        }
    }
    std::remove(file.c_str());
    std::remove(dir.c_str());
}

TEST_CASE(
    "NucleusSampler::readConfigurationFile picks the deuteron file matching "
    "the requested Jz, or a random one for an unpolarized deuteron") {
    Parameters param;
    makeInitTestParam(param, 4);
    const std::string dir = "test_deuteron_configs_tmp";
    REQUIRE(std::system(("mkdir -p " + dir).c_str()) == 0);
    // one configuration of 2 nucleons each; all entries tag the file
    const std::string pol0 = dir + "/DeuteronPol0Configs.bin.in";
    const std::string polpm1 = dir + "/DeuteronPolpm1Configs.bin.in";
    for (const auto &[file, tag] :
         {std::pair<std::string, float> {pol0, 0.f}, {polpm1, 1.f}}) {
        std::ofstream out(file, std::ios::binary);
        for (int k = 0; k < 2 * 3; k++) {
            out.write(reinterpret_cast<const char *>(&tag), sizeof(tag));
        }
    }
    param.nucleus.nuclearConfigurationsPath = dir;

    for (const auto &[Jz, expectedTag] :
         {std::pair<double, float> {0., 0.f}, {1., 1.f}, {-1., 1.f}}) {
        CAPTURE(Jz);
        NucleusSampler sampler;
        // polarizationFlag != 0: the file is chosen by Jz, not randomly
        const std::vector<std::vector<float>> configs =
            sampler.readConfigurationFile(2, 0, 1, Jz, nullptr, &param);
        REQUIRE(configs.size() == 1);
        CHECK(configs[0][0] == expectedTag);
    }

    // polarizationFlag 0: Jz = 0 with probability 1/3, else |Jz| = 1,
    // whatever Jz is given
    Random random;
    random.init_genrand64(23ULL);
    const int samples = 300;
    int pol0Picks = 0;
    for (int i = 0; i < samples; i++) {
        NucleusSampler sampler;
        const std::vector<std::vector<float>> configs =
            sampler.readConfigurationFile(2, 0, 0, 1., &random, &param);
        REQUIRE(configs.size() == 1);
        if (configs[0][0] == 0.f) pol0Picks++;
    }
    CHECK(
        static_cast<double>(pol0Picks) / samples
        == doctest::Approx(1. / 3.).epsilon(0.25));
    std::remove(pol0.c_str());
    std::remove(polpm1.c_str());
    std::remove(dir.c_str());
}

TEST_CASE(
    "configurationFileName and lightNucleusOptionError: the files by mass "
    "number and lightNucleusOption") {
    // species with several files
    CHECK(configurationFileName(3, 1) == "triton.bin.in");
    CHECK(
        configurationFileName(16, 4)
        == "O16_NLEFT_dmin0.5fm_positiveweights.bin.in");
    CHECK(configurationFileName(20, 0) == configurationFileName(20, 2));
    CHECK(configurationFileName(40, 4) == "Ar40_NLEFT.bin.in");
    for (const auto &[nucleusA, option] :
         std::vector<std::pair<int, int>> {{3, 0}, {16, 5}, {20, 3}, {40, 0}}) {
        CAPTURE(nucleusA);
        CAPTURE(option);
        CHECK(lightNucleusOptionError(nucleusA, option).empty());
    }
    // options without a file are rejected
    CHECK(configurationFileName(40, 5).empty());
    CHECK(configurationFileName(20, 1).empty());
    const std::string error = lightNucleusOptionError(40, 5);
    CHECK(error.find("A = 40") != std::string::npos);
    CHECK(error.find("allowed values: 0 4") != std::string::npos);
    CHECK_FALSE(lightNucleusOptionError(12, 2).empty());
    // species with one file, or none, ignore the option
    CHECK(configurationFileName(208, 5) == "Pb208.bin.in");
    CHECK(lightNucleusOptionError(208, 5).empty());
    CHECK(lightNucleusOptionError(1, 3).empty());
    CHECK(configurationFileName(63, 0).empty());
}
