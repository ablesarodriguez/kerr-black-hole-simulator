#ifndef CAMERA_H
#define CAMERA_H

struct Camera {
    // Camera position in Schwarzschild coordinates: [t, r, theta, phi]
    double position[4];
    
    // Viewport resolution and lens settings
    int screen_width;
    int screen_height;
    double field_of_view; // Field of view in radians
    
    // Constructor to setup the observer
    Camera(int width, int height, double fov, double black_hole_mass) {
        screen_width = width;
        screen_height = height;
        field_of_view = fov;
        
        // Initialize camera coordinates far away from the black hole
        
        // position[0] -> t = 0.0 (Initial coordinate time)
        position[0] = 0.0;
        
        // position[1] -> r = 100M (Far enough to approximate flat spacetime)
        position[1] = 100.0 * black_hole_mass;
        
        // position[2] -> theta = pi/2 (Placed exactly on the equatorial plane)
        position[2] = 3.14159265359 / 2.0;
        
        // position[3] -> phi = 0.0 (Starting azimuthal angle)
        position[3] = 0.0;
    }
};

#endif // CAMERA_H