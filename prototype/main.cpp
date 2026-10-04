#include <iostream>
#include <fstream>
#include <cmath>

#include "BlackHole.h"
#include "Vector3.h"

int main() {
    double mass = 1.0;
    BlackHole bh(mass);
    
    // Increase resolution for the 3D render
    int width = 800;
    int height = 800;
    
    // 3D Camera Setup
    // Positioned slightly above the equator (y=2.0) to see the disk in perspective
    Vector3 cam_pos(0.0, 2.0, 20.0); 
    
    std::ofstream image_file("black_hole_3d.ppm");
    image_file << "P3\n" << width << " " << height << "\n255\n";

    std::cout << "Starting 3D relativistic ray marching..." << std::endl;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            
            // Map screen to 3D viewing plane
            double normalized_x = (2.0 * x / (double)width) - 1.0;
            double normalized_y = 1.0 - (2.0 * y / (double)height);
            
            // Set initial ray direction into the screen (-Z)
            Vector3 ray_dir = Vector3(normalized_x, normalized_y, -1.0).normalize();
            
            Vector3 pos = cam_pos;
            Vector3 vel = ray_dir;
            
            // Ray Marching parameters
            double step_size = 0.1;
            int max_steps = 1000;
            
            // State flags
            bool hit_black_hole = false;
            bool hit_disk = false;
            double disk_radius_hit = 0.0;

            // Integration Loop (Moving the photon step by step)
            for (int steps = 0; steps < max_steps; steps++) {
                double r = pos.length();
                
                // 1. Check if photon fell into the Event Horizon
                if (r < bh.r_horizon) {
                    hit_black_hole = true;
                    break;
                }
                
                // 2. Check if photon escaped far away
                if (r > 30.0) {
                    break; 
                }
                
                // 3. 3D Accretion Disk Check (The plane where Y = 0)
                // If the Y coordinate signs are different, we crossed the flat disk
                double next_y = pos.y + vel.y * step_size;
                if (pos.y * next_y <= 0.0) {
                    // Check if the crossing happened within the disk's physical limits (3M to 10M)
                    if (r > bh.r_photon_sphere && r < 10.0 * bh.M) {
                        hit_disk = true;
                        disk_radius_hit = r;
                        break; // Stop ray, we see the disk!
                    }
                }
                
                // 4. Relativistic Gravity Calculation (The magic of General Relativity)
                // Acceleration = - (3M * h^2 / r^5) * position
                Vector3 h = pos.cross(vel); // Angular momentum vector
                double h2 = h.dot(h);
                Vector3 accel = pos * (-3.0 * bh.M * h2 / std::pow(r, 5));
                
                // Update velocity and position
                vel = vel + accel * step_size;
                vel = vel.normalize(); // Photons always travel at speed of light (length 1)
                pos = pos + vel * step_size;
            }

            // Coloring the pixel based on what the photon hit
            int r_col = 0, g_col = 0, b_col = 0;

            if (hit_black_hole) {
                // Absolute darkness
                r_col = 0; g_col = 0; b_col = 0;
            } else if (hit_disk) {
                // The Accretion Disk - Hotter (brighter) closer to the center
                double intensity = 1.0 - ((disk_radius_hit - bh.r_photon_sphere) / (10.0 * bh.M - bh.r_photon_sphere));
                
                r_col = 255;
                g_col = (int)(150 * intensity) + 50;
                b_col = (int)(50 * intensity);
            } else {
                // Deep space background
                r_col = 5; g_col = 10; b_col = 20;
            }

            image_file << r_col << " " << g_col << " " << b_col << "\n";
        }
    }

    image_file.close();
    std::cout << "3D Render complete! Open 'black_hole_3d.ppm'" << std::endl;

    return 0;
}