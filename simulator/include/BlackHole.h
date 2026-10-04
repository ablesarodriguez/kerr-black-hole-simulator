#ifndef BLACKHOLE_H
#define BLACKHOLE_H

struct BlackHole {
    double M;                // Mass of the black hole
    double r_horizon;        // Event horizon (2M)
    double r_photon_sphere;  // Photon sphere (3M)
    double r_isco;           // Innermost stable circular orbit (6M)

    BlackHole(double mass) {
        M = mass;
        r_horizon = 2.0 * M;
        r_photon_sphere = 3.0 * M;
        r_isco = 6.0 * M;
    }
};

#endif // BLACKHOLE_H