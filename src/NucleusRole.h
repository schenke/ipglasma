// NucleusRole.h is part of the IP-Glasma solver.

#ifndef SRC_NUCLEUSROLE_H_
#define SRC_NUCLEUSROLE_H_

/**
 * Which of the two colliding nuclei a lattice/nucleon-sampling operation
 * applies to.
 */
enum class NucleusRole {
    /// The first (\f$A\f$) nucleus.
    Projectile,
    /// The second (\f$B\f$) nucleus.
    Target,
};

#endif  // SRC_NUCLEUSROLE_H_
