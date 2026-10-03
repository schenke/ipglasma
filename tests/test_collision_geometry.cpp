#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "CollisionGeometry.h"
#include "Glauber.h"
#include "Lattice.h"
#include "Parameters.h"
#include "PhysConst.h"
#include "Random.h"
#include "doctest.h"
#include "test_helpers.h"

namespace {
ReturnValue nucleonAt(double x, double y, bool proton = true) {
    ReturnValue n {};
    n.x = x;
    n.y = y;
    n.proton = proton;
    return n;
}

void removeGeometryFiles() {
    std::remove("NcollList0.dat");
    std::remove("NpartList0.dat");
    std::remove("NgluonEstimators0.dat");
    std::remove("usedParameters0.dat");
}

int countLines(const std::string &fileName) {
    std::ifstream in(fileName);
    int n = 0;
    for (std::string line; std::getline(in, line);) {
        if (!line.empty()) n++;
    }
    return n;
}
}  // namespace

TEST_CASE(
    "CollisionGeometry::computeNcollList (hard-sphere mode) marks a p+p "
    "pair collided and counts it once") {
    Parameters param;
    param.collision.sigmaNN = 4.2;
    param.collision.gaussianWounding = false;
    param.event.eventId = 0;
    std::vector<ReturnValue> nucleusA {nucleonAt(0., 0.)};
    std::vector<ReturnValue> nucleusB {nucleonAt(0., 0.)};
    CollisionGeometry geometry(nucleusA, nucleusB);

    const double d2 = param.collision.sigmaNN / (M_PI * 10.);  // in fm^2
    int Ncoll = 0;
    geometry.computeNcollList(
        &param, nullptr, d2, /*b=*/0.0, /*phiRP=*/0.0, Ncoll);

    CHECK(Ncoll == 1);
    CHECK(nucleusA[0].collided == 1);
    CHECK(nucleusB[0].collided == 1);
    removeGeometryFiles();
}

TEST_CASE(
    "CollisionGeometry::determineNpartAndNcoll counts the wounded nucleons "
    "of both nuclei at the sampled impact parameter") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.collision.sigmaNN = 42.;  // d = sqrt(4.2 / pi) = 1.16 fm
    param.collision.gaussianWounding = false;
    // A sits at +b/2, B at -b/2 along x: A's nucleon at x = -1 meets B's
    // two nucleons at x = 1 and x = 1.5, A's nucleon at x = 5 meets nobody
    param.event.b = 2.;
    param.event.phiRP = 0.;
    std::vector<ReturnValue> nucleusA {nucleonAt(-1., 0.), nucleonAt(5., 0.)};
    std::vector<ReturnValue> nucleusB {
        nucleonAt(1., 0.), nucleonAt(1.5, 0.3, false), nucleonAt(-6., 0.)};
    CollisionGeometry geometry(nucleusA, nucleusB);

    int Npart = 0, Ncoll = 0;
    REQUIRE(geometry.determineNpartAndNcoll(&param, nullptr, Npart, Ncoll));
    CHECK(Ncoll == 2);
    CHECK(Npart == 3);
    CHECK(param.event.Npart == 3);
    CHECK(nucleusA[0].collided == 1);
    CHECK(nucleusA[1].collided == 0);
    CHECK(nucleusB[0].collided == 1);
    CHECK(nucleusB[1].collided == 1);
    CHECK(nucleusB[2].collided == 0);
    CHECK(countLines("NcollList0.dat") == 2);
    CHECK(countLines("NpartList0.dat") == 5);

    // a fixed Npart that doesn't match rejects the event
    param.collision.useFixedNpart = 4;
    param.event.success = 1;
    CHECK_FALSE(geometry.determineNpartAndNcoll(&param, nullptr, Npart, Ncoll));
    CHECK(param.event.success == 0);
    removeGeometryFiles();
}

TEST_CASE(
    "CollisionGeometry::determineNpartAndNcoll treats a p+p event as "
    "collided even without overlap") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.collision.sigmaNN = 42.;
    param.event.b = 10.;
    std::vector<ReturnValue> nucleusA {nucleonAt(0., 0.)};
    std::vector<ReturnValue> nucleusB {nucleonAt(0., 0.)};
    CollisionGeometry geometry(nucleusA, nucleusB);
    int Npart = 0, Ncoll = 0;
    REQUIRE(geometry.determineNpartAndNcoll(&param, nullptr, Npart, Ncoll));
    CHECK(Ncoll == 0);
    CHECK(Npart == 2);
    removeGeometryFiles();
}

TEST_CASE(
    "CollisionGeometry::sampleImpactParameter samples b in [bMin, bMax] and "
    "resets the collided flags") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.collision.useNucleus = true;
    param.collision.bMin = 2.;
    param.collision.bMax = 3.;
    param.collision.rotateReactionPlane = false;
    std::vector<ReturnValue> nucleusA {nucleonAt(0., 0.)};
    std::vector<ReturnValue> nucleusB {nucleonAt(0., 0.), nucleonAt(1., 0.)};
    nucleusA[0].collided = 1;
    nucleusB[1].collided = 1;
    CollisionGeometry geometry(nucleusA, nucleusB);
    Random random;
    random.init_genrand64(4);

    for (bool linear : {true, false}) {
        CAPTURE(linear);
        param.collision.sampleBFromLinearDistribution = linear;
        for (int i = 0; i < 20; i++) {
            geometry.sampleImpactParameter(&param, &random);
            CHECK(param.event.b >= 2.);
            CHECK(param.event.b <= 3.);
            CHECK(param.event.phiRP == 0.);
        }
    }
    CHECK(nucleusA[0].collided == 0);
    CHECK(nucleusB[1].collided == 0);

    param.collision.rotateReactionPlane = true;
    geometry.sampleImpactParameter(&param, &random);
    CHECK(param.event.phiRP >= 0.);
    CHECK(param.event.phiRP < 2. * M_PI);

    param.collision.useNucleus = false;
    geometry.sampleImpactParameter(&param, &random);
    CHECK(param.event.b == 0.);
}

TEST_CASE(
    "CollisionGeometry::computeQuantities rejects an event without overlap "
    "region") {
    const int N = 8;
    Parameters param;
    makeInitTestParam(param, N);
    param.collision.sigmaNN = 42.;
    param.event.b = 0.;
    param.event.success = 1;
    // nucleons far apart: the p+p rule wounds both, but no lattice cell is
    // within the wounding distance of both
    std::vector<ReturnValue> nucleusA {nucleonAt(-3., 0.)};
    std::vector<ReturnValue> nucleusB {nucleonAt(3., 0.)};
    CollisionGeometry geometry(nucleusA, nucleusB);
    Lattice lat(&param, N);
    std::remove("NgluonEstimators0.dat");
    geometry.computeQuantities(&lat, &param, nullptr);
    CHECK(param.event.success == 0);
    CHECK(param.event.averageQs == 0.);
    // the no-overlap branch returns before the estimators are written
    CHECK_FALSE(std::ifstream("NgluonEstimators0.dat").good());
    removeGeometryFiles();
}

TEST_CASE(
    "CollisionGeometry::computeAndSetRunningAlphaS respects nFlavors/LambdaQCD "
    "instead "
    "of assuming 3 flavors and LambdaQCD=0.2, and uses the same muZero/c "
    "as computeRunningCouplingGfactor()") {
    Parameters param;
    makeInitTestParam(param, 4);
    param.coupling.runningCoupling = true;
    param.coupling.runWithKt = false;
    param.coupling.runWithQs = 1;  // average Qs
    param.coupling.runningCouplingQsFactor = 0.5;
    param.event.averageQsAvg = 1.3;
    param.coupling.mu0 = 0.3;
    param.coupling.c = 0.2;

    std::vector<ReturnValue> nucleusA, nucleusB;
    CollisionGeometry geometry(nucleusA, nucleusB);

    param.coupling.nFlavors = 3;
    param.coupling.LambdaQCD = 0.2;
    geometry.computeAndSetRunningAlphaS(&param);
    CHECK(param.event.alphas == doctest::Approx(0.5922901359561648));

    param.coupling.nFlavors = 4;
    param.coupling.LambdaQCD = 0.25;
    geometry.computeAndSetRunningAlphaS(&param);
    CHECK(param.event.alphas == doctest::Approx(0.789051392062008));
}

TEST_CASE(
    "CollisionGeometry::scanOverlap reads the shifted nuclei only on the "
    "lattice, including the corner cell (0, 0)") {
    const int N = 6;
    Parameters param;
    makeInitTestParam(param, N);  // a = 1 fm
    param.coupling.g = 1.;
    param.colorCharge.QsMuRatio = 1.;
    param.nucleus.useSmoothNucleus = false;
    Lattice lat(&param, N);
    for (int pos = 0; pos < N * N; pos++) {
        lat.cells[pos]->setg2mu2A(1.);
        lat.cells[pos]->setg2mu2B(1.);
    }
    std::vector<ReturnValue> nucleusA, nucleusB;
    CollisionGeometry geometry(nucleusA, nucleusB);
    const double a = 1.;
    // Q_s^2 of one cell with g^2 mu^2 = 1, as summed into averageQs2min2
    const double cellQs2 = PhysConst::hbarc * PhysConst::hbarc;

    auto fullLatticeSum = [&](double b, double phiRP) {
        double averageQs = 0., averageQs2 = 0., averageQs2Avg = 0.,
               averageQs2min = 0., averageQs2min2 = 0., Tpp = 0.;
        int count = 0;
        geometry.scanOverlap(
            &lat, &param, N, a, b, phiRP, averageQs, averageQs2, averageQs2Avg,
            averageQs2min, averageQs2min2, Tpp, count);
        return averageQs2min2 / cellQs2;
    };
    // without a shift every cell, including (0, 0), has both nuclei
    CHECK(fullLatticeSum(0., 0.) == doctest::Approx(N * N));
    // shifting by b = 2 fm along y moves A by -1 and B by +1 cell in y; a
    // cell only counts if both shifted cells lie on the lattice, which
    // leaves N - 2 rows (a flattened index used to let (ix, -1) alias
    // (ix - 1, N - 1))
    CHECK(fullLatticeSum(2., M_PI / 2.) == doctest::Approx(N * (N - 2)));
}
