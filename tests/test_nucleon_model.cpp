#include <algorithm>
#include <array>
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

TEST_CASE(
    "HotSpotNucleon: NqFluc adds a Poisson fluctuation to the number of hot "
    "spots") {
    Parameters param;
    makeHotSpotParam(param);
    param.subnucleon.NqFluc = 0.5;
    const HotSpotNucleon model(param);
    Random random;
    random.init_genrand64(13ULL);
    random.gslRandomInit(13ULL);

    const int samples = 20000;
    double sum = 0., sumSq = 0.;
    for (int i = 0; i < samples; ++i) {
        const int n = model.sampleNumberOfPartons(random);
        CHECK(n >= 3);
        sum += n;
        sumSq += n * n;
    }
    const double mean = sum / samples;
    CHECK(mean == doctest::Approx(3.5).epsilon(0.02));
    CHECK(sumSq / samples - mean * mean == doctest::Approx(0.5).epsilon(0.05));
}

TEST_CASE(
    "HotSpotNucleon: shiftConstituentQuarkProtonOrigin puts the hot spots' "
    "center of mass at the nucleon center") {
    for (bool shift : {true, false}) {
        CAPTURE(shift);
        Parameters param;
        makeHotSpotParam(param);
        param.subnucleon.shiftConstituentQuarkProtonOrigin = shift;
        const HotSpotNucleon model(param);
        Random random;
        random.init_genrand64(17ULL);

        double largestOffset = 0.;
        for (int i = 0; i < 20; ++i) {
            const std::unique_ptr<NucleonProfile> profile =
                model.sample(random, nucleonAt(0.4, 0.2));
            const auto *hotSpots =
                dynamic_cast<const HotSpotProfile *>(profile.get());
            REQUIRE(hotSpots != nullptr);
            double meanX = 0., meanY = 0.;
            for (int iq = 0; iq < hotSpots->numberOfHotSpots(); ++iq) {
                meanX += hotSpots->offsetsX()[iq];
                meanY += hotSpots->offsetsY()[iq];
            }
            meanX /= hotSpots->numberOfHotSpots();
            meanY /= hotSpots->numberOfHotSpots();
            largestOffset = std::max(largestOffset, std::hypot(meanX, meanY));
        }
        if (shift) {
            CHECK(largestOffset < 1e-12);
        } else {
            CHECK(largestOffset > 0.05);  // fm
        }
    }
}

TEST_CASE(
    "HotSpotNucleon: dqMin makes hot spots closer than dqMin rare (omega != 1 "
    "places them in the transverse plane)") {
    // The minimum distance is best effort: a hot spot keeps its sampled
    // radius and only its angle is redrawn (up to 100 times), so two hot
    // spots at nearly the same small radius can stay closer than dqMin.
    auto closePairFraction = [](double dqMin) {
        Parameters param;
        makeHotSpotParam(param);
        param.subnucleon.omega = 2.;
        param.subnucleon.dqMin = dqMin;
        const HotSpotNucleon model(param);
        Random random;
        random.init_genrand64(19ULL);
        random.setGammaIncCDF(param.subnucleon.omega);
        int pairs = 0, close = 0;
        for (int i = 0; i < 2000; ++i) {
            const std::unique_ptr<NucleonProfile> profile =
                model.sample(random, nucleonAt(0., 0.));
            const auto *hotSpots =
                dynamic_cast<const HotSpotProfile *>(profile.get());
            REQUIRE(hotSpots != nullptr);
            const std::vector<double> &x = hotSpots->offsetsX();
            const std::vector<double> &y = hotSpots->offsetsY();
            for (std::size_t a = 0; a < x.size(); ++a) {
                for (std::size_t b = 0; b < a; ++b) {
                    pairs++;
                    if (std::hypot(x[a] - x[b], y[a] - y[b]) < 0.3) close++;
                }
            }
        }
        return static_cast<double>(close) / pairs;
    };
    const double without = closePairFraction(0.);
    const double with = closePairFraction(0.3);
    CHECK(without > 0.1);
    CHECK(with < without / 4.);
}

TEST_CASE(
    "HotSpotNucleon: reference thicknesses with every option on (omega, "
    "dqMin, NqFluc, BGqVar, smearQs)") {
    // Values of the current code, whose full runs are byte-identical to
    // the code before the nucleon-model interface; they change if the
    // sampling (including the order of the random draws) changes.
    const double expected[2][3][3] = {
        // omega = 1 (3D Gaussian positions)
        {{0.06130905822699826, 0.0030647365940887714, 0.00013074215615949835},
         {0.0012159547566270032, 0.1424820795261299, 0.0059570302662717454},
         {0.012469481836384959, 0.095459167771926923, 0.044697293222758824}},
        // omega = 2 (transverse gamma-distributed radii)
        {{0.14102867950110543, 0.0031792272509544805, 0.001158864499172682},
         {0.031111785576826728, 5.8065139259764331e-06, 2.2767536952053522e-09},
         {0.049810738087300084, 0.011669559549845221, 0.0016938040070285157}}};
    const int expectedNq[2][3] = {{3, 3, 3}, {3, 2, 4}};
    const double points[3][2] = {{0.3, -0.2}, {0.6, 0.1}, {-0.1, -0.5}};

    for (int io = 0; io < 2; ++io) {
        const double omega = io == 0 ? 1. : 2.;
        CAPTURE(omega);
        Parameters param;
        makeHotSpotParam(param);
        param.subnucleon.BGqVar = 0.1;
        param.subnucleon.dqMin = 0.2;
        param.subnucleon.omega = omega;
        param.subnucleon.NqBase = 2.5;
        param.subnucleon.NqFluc = 0.5;
        param.subnucleon.smearQs = true;
        param.subnucleon.smearingWidth = 0.5;
        Random random;
        random.init_genrand64(2024ULL);
        random.gslRandomInit(2024ULL);
        random.setGammaIncCDF(omega);
        const HotSpotNucleon model(param);
        for (int i = 0; i < 3; ++i) {
            CAPTURE(i);
            const std::unique_ptr<NucleonProfile> profile =
                model.sample(random, nucleonAt(0.3, -0.2));
            const auto *hotSpots =
                dynamic_cast<const HotSpotProfile *>(profile.get());
            REQUIRE(hotSpots != nullptr);
            CHECK(hotSpots->numberOfHotSpots() == expectedNq[io][i]);
            for (int ip = 0; ip < 3; ++ip) {
                CHECK(
                    profile->thickness(points[ip][0], points[ip][1])
                    == doctest::Approx(expected[io][i][ip]).epsilon(1e-12));
            }
        }
    }
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
    param.subnucleon.nucleonModel = "strings";
    CHECK(
        dynamic_cast<StringyNucleon *>(NucleonModel::create(param).get())
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

namespace {
double distance(const Point3 &a, const Point3 &b) {
    return std::sqrt(
        (a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1])
        + (a[2] - b[2]) * (a[2] - b[2]));
}

double totalLength(const Point3 &p, const std::array<Point3, 3> &points) {
    return distance(p, points[0]) + distance(p, points[1])
           + distance(p, points[2]);
}

// Whether p is (to rounding) where the total length to the points is
// smallest: no step of size h in any of the 26 directions shortens it.
bool isMinimum(const Point3 &p, const std::array<Point3, 3> &points) {
    const double h = 1e-6;
    const double length = totalLength(p, points);
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                const Point3 q {p[0] + h * dx, p[1] + h * dy, p[2] + h * dz};
                if (totalLength(q, points) < length - 1e-12) return false;
            }
        }
    }
    return true;
}
}  // namespace

TEST_CASE(
    "StringyNucleon::fermatPoint: the centroid of an equilateral triangle, "
    "the vertex of an angle of at least 120 degrees, a coincident point") {
    const double s3 = std::sqrt(3.);
    const std::array<Point3, 3> equilateral {
        Point3 {1., 0., 0.5}, Point3 {-0.5, s3 / 2., 0.5},
        Point3 {-0.5, -s3 / 2., 0.5}};
    const Point3 centroid = StringyNucleon::fermatPoint(equilateral);
    CHECK(distance(centroid, Point3 {0., 0., 0.5}) < 1e-12);

    // an angle of 150 degrees at the first point
    const double c = std::cos(150. * M_PI / 180.);
    const double s = std::sin(150. * M_PI / 180.);
    const std::array<Point3, 3> obtuse {
        Point3 {0.1, 0.2, 0.3}, Point3 {1.1, 0.2, 0.3},
        Point3 {0.1 + 0.7 * c, 0.2 + 0.7 * s, 0.3}};
    CHECK(StringyNucleon::fermatPoint(obtuse) == obtuse[0]);
    // exactly 120 degrees: the vertex is also where the sum is smallest
    const std::array<Point3, 3> limit {
        Point3 {0., 0., 0.}, Point3 {1., 0., 0.}, Point3 {-0.5, s3 / 2., 0.}};
    CHECK(distance(StringyNucleon::fermatPoint(limit), limit[0]) < 1e-12);

    const std::array<Point3, 3> coincident {
        Point3 {0.3, 0.1, 0.}, Point3 {-0.2, 0.4, 0.1}, Point3 {0.3, 0.1, 0.}};
    CHECK(StringyNucleon::fermatPoint(coincident) == coincident[0]);
}

TEST_CASE(
    "StringyNucleon::fermatPoint minimizes the total distance for random "
    "triangles in 3D, independently of the order of the points") {
    Random random;
    random.init_genrand64(29ULL);
    int vertexCases = 0, interiorCases = 0;
    for (int i = 0; i < 500; ++i) {
        std::array<Point3, 3> points;
        for (Point3 &p : points) {
            for (double &coordinate : p) coordinate = random.gauss(0., 0.5);
        }
        CAPTURE(i);
        const Point3 fermat = StringyNucleon::fermatPoint(points);
        CHECK(isMinimum(fermat, points));
        bool atVertex = false;
        for (const Point3 &p : points) {
            if (distance(fermat, p) < 1e-12) atVertex = true;
        }
        if (atVertex) {
            vertexCases++;
        } else {
            interiorCases++;
            // each pair of points is seen at 120 degrees
            for (int a = 0; a < 3; a++) {
                const int b = (a + 1) % 3;
                double dotProduct = 0.;
                for (int k = 0; k < 3; k++) {
                    dotProduct +=
                        (points[a][k] - fermat[k]) * (points[b][k] - fermat[k]);
                }
                CHECK(
                    dotProduct
                        / (distance(points[a], fermat)
                           * distance(points[b], fermat))
                    == doctest::Approx(-0.5).epsilon(1e-9));
            }
        }
        const std::array<Point3, 3> permuted {points[2], points[0], points[1]};
        CHECK(distance(StringyNucleon::fermatPoint(permuted), fermat) < 1e-12);
    }
    // both cases occur for random triangles
    CHECK(vertexCases > 10);
    CHECK(interiorCases > 10);
}

TEST_CASE(
    "StringyNucleon::sampleStrings moves each of the three hot spots to a "
    "point on its string from the junction") {
    for (double omega : {1., 2.}) {
        CAPTURE(omega);
        Parameters param;
        makeHotSpotParam(param);
        param.subnucleon.nucleonModel = "strings";
        param.subnucleon.omega = omega;
        param.subnucleon.NqBase = 5.;  // not used: always three hot spots
        const StringyNucleon model(param);
        Random random;
        random.init_genrand64(31ULL);
        random.setGammaIncCDF(omega);

        bool someZ = false;
        for (int n = 0; n < 50; ++n) {
            const StringConfiguration strings = model.sampleStrings(random);
            const HotSpotConfiguration &hotSpots = strings.hotSpots;
            REQUIRE(hotSpots.x.size() == 3);
            REQUIRE(strings.t.size() == 3);
            REQUIRE(strings.x.size() == 3);
            std::array<Point3, 3> ends;
            for (int i = 0; i < 3; i++) {
                ends[i] = {hotSpots.x[i], hotSpots.y[i], hotSpots.z[i]};
                if (std::abs(hotSpots.z[i]) > 1e-3) someZ = true;
            }
            const Point3 &J = strings.junction;
            CHECK(distance(J, StringyNucleon::fermatPoint(ends)) < 1e-15);
            for (int i = 0; i < 3; i++) {
                const double t = strings.t[i];
                CHECK(t > 0.);
                CHECK(t < 1.);
                CHECK(
                    strings.x[i]
                    == doctest::Approx(J[0] + t * (ends[i][0] - J[0])));
                CHECK(
                    strings.y[i]
                    == doctest::Approx(J[1] + t * (ends[i][1] - J[1])));
            }
        }
        // omega 1 samples the hot spots in 3D, omega != 1 in the plane
        CHECK(someZ == (omega == 1.));
    }
}

TEST_CASE(
    "StringyNucleon: three moved hot spots whose profile integrates to the "
    "Qs normalization") {
    Parameters param;
    makeHotSpotParam(param);
    param.subnucleon.nucleonModel = "strings";
    const StringyNucleon model(param);
    Random random;
    random.init_genrand64(37ULL);
    for (int i = 0; i < 5; ++i) {
        const std::unique_ptr<NucleonProfile> profile =
            model.sample(random, nucleonAt(-0.4, 0.7));
        const auto *hotSpots =
            dynamic_cast<const HotSpotProfile *>(profile.get());
        REQUIRE(hotSpots != nullptr);
        CHECK(hotSpots->numberOfHotSpots() == 3);
        CHECK(
            integrate(*profile, -0.4, 0.7)
            == doctest::Approx(1.).epsilon(1e-3));
    }
}

TEST_CASE(
    "StringyNucleon: reference thicknesses with every option on (omega, "
    "dqMin, BGqVar, smearQs)") {
    // Values of the current code (its geometry is checked by the tests
    // above); they change if the sampling, including the order of the
    // random draws, changes.
    const double expected[2][3][3] = {
        // omega = 1 (hot spots and junction in 3D)
        {{0.07117971749068952, 0.0097705399534210192, 0.00073540124885799021},
         {0.067710670547497645, 0.0018023507228984231, 8.3750080501346203e-06},
         {0.04893515017337597, 0.006223252893511509, 0.0072116750964462545}},
        // omega = 2 (in the transverse plane)
        {{0.16727568843466295, 0.016035469215763525, 0.0010477318473265445},
         {0.023444465761166905, 0.0042044758566608734, 0.0016717585390443976},
         {0.014312410429116873, 5.5722266360508553e-06,
          2.7781712707968981e-13}}};
    const double points[3][2] = {{0.3, -0.2}, {0.6, 0.1}, {-0.1, -0.5}};

    for (int io = 0; io < 2; ++io) {
        const double omega = io == 0 ? 1. : 2.;
        CAPTURE(omega);
        Parameters param;
        makeHotSpotParam(param);
        param.subnucleon.nucleonModel = "strings";
        param.subnucleon.BGqVar = 0.1;
        param.subnucleon.dqMin = 0.2;
        param.subnucleon.omega = omega;
        param.subnucleon.smearQs = true;
        param.subnucleon.smearingWidth = 0.5;
        Random random;
        random.init_genrand64(2024ULL);
        random.gslRandomInit(2024ULL);
        random.setGammaIncCDF(omega);
        const StringyNucleon model(param);
        for (int i = 0; i < 3; ++i) {
            CAPTURE(i);
            const std::unique_ptr<NucleonProfile> profile =
                model.sample(random, nucleonAt(0.3, -0.2));
            for (int ip = 0; ip < 3; ++ip) {
                CHECK(
                    profile->thickness(points[ip][0], points[ip][1])
                    == doctest::Approx(expected[io][i][ip]).epsilon(1e-12));
            }
        }
    }
}
