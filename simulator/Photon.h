#ifndef PHOTON_H
#define PHOTON_H

struct Photon {
    // Schwarzschild coordinates: [t, r, theta, phi]
    double pos[4];
    
    // Four-momentum: rate of change for each coordinate
    double p[4];

    // Constants of motion
    double e;  // Energy 
    double l;  // Angular momentum 
    double b;  // Impact parameter (b = l / e)
    
    // State flag for the rendering loop
    bool is_alive;  
};

#endif // PHOTON_H