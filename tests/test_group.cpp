#include "doctest.h"

#include "Group.h"

TEST_CASE("Group: fundamental-rep generators are Hermitian and traceless") {
    Group group;
    for (int a = 0; a < 8; ++a) {
        Matrix t = group.getT(a);
        Matrix adjoint = t;
        adjoint.conjg();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                CHECK(std::abs(t.get(i, j) - adjoint.get(i, j)) < 1e-14);
            }
        }
        CHECK(std::abs(t.trace()) < 1e-14);
    }
}

TEST_CASE("Group: const getT() overload returns the same generators") {
    Group group;
    const Group &constGroup = group;
    for (int a = 0; a < 8; ++a) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                CHECK(
                    std::abs(
                        constGroup.getT(a).get(i, j) - group.getT(a).get(i, j))
                    < 1e-14);
            }
        }
    }
}

TEST_CASE("Group: generators are normalized to Tr(t^a t^b) = 0.5 delta^ab") {
    Group group;
    for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
            Matrix product = group.getT(a) * group.getT(b);
            const double expected = (a == b) ? 0.5 : 0.0;
            CHECK(std::abs(product.trace() - expected) < 1e-13);
        }
    }
}
