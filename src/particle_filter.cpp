/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "helper_functions.h"
#include "particle_filter.h"

using std::string;
using std::vector;
using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  if (is_initialized){
	return;
  }

  num_particles = 100;  // TODO: Set the number of particles
  particles.resize(num_particles);
  weights.resize(num_particles);

  std::default_random_engine gen;

  // Creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  for (int i = 0; i < num_particles; i++){
    Particle par;
    par.id = i;
    par.x = dist_x(gen);
    par.y = dist_y(gen);
    par.theta = dist_theta(gen);
    par.weight = 1.0;
   
    particles.push_back(par);
    weights[i] = 1.0;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
   
   std::default_random_engine gen;
   normal_distribution<double> dist_x(0, std_pos[0]);
   normal_distribution<double> dist_y(0, std_pos[1]);
   normal_distribution<double> dist_theta(0, std_pos[2]);

   // In case of unchanging yaw rate (division by zero will happen)
   for (int i = 0; i < num_particles; i++){
     if (fabs(yaw_rate) < 0.00001){
        particles[i].x += velocity * delta_t * cos(particles[i].theta);
        particles[i].y += velocity * delta_t * sin(particles[i].theta);
     }
     else {
       particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
       particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
       particles[i].theta += yaw_rate * delta_t;
     }

     // Adding noise
     particles[i].x += dist_x(gen);
     particles[i].y += dist_y(gen);
     particles[i].theta += dist_theta(gen);
   }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
   
   // filtering out the landmarks which are out of range from particles view
  
   double x_p, y_p, theta_p, weight, x_m, y_m, x_c, y_c;
   double mu_x, mu_y, expon;
   double sig_x = std_landmark[0];
   double sig_y = std_landmark[1];
   double gauss_norm = 1 / (2*M_PI*sig_x*sig_y);

   for (int i = 0; i < num_particles; i++){

     x_p = particles[i].x;
     y_p = particles[i].y;
     theta_p = particles[i].theta;
     weight = 1.0;

     vector<LandmarkObs> nearby_lms;
     LandmarkObs nearby_lm;

     //Getting the landmarks in sensor range
     for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++){
       double lm_x = map_landmarks.landmark_list[j].x_f;
       double lm_y = map_landmarks.landmark_list[j].y_f;
       int lm_id = map_landmarks.landmark_list[j].id_i;
       double lm_dist = dist(x_p, y_p, lm_x, lm_y);
       
       if (lm_dist < sensor_range){
         nearby_lm.id = lm_id;
         nearby_lm.x = lm_x;
         nearby_lm.y = lm_y;
         nearby_lms.push_back(nearby_lm);
       }		
     }

     // Transform from vehicle to particle coordinates.
     //vector<LandmarkObs> transformed_obs;
     for(unsigned int k = 0; k < observations.size(); k++) {
       x_c = observations[k].x;
       y_c = observations[k].y;
       x_m = cos(theta_p) * x_c - sin(theta_p) * y_c + x_p;
       y_m = sin(theta_p) * x_c + cos(theta_p) * y_c + y_p;

       double obs_min_dist = 999;
       int min_obs_index = 0;
       
       for (unsigned int l = 0; l < nearby_lms.size(); l++){
         double obs_to_lm = dist(x_m, y_m, nearby_lms[l].x, nearby_lms[l].y);
         if (obs_to_lm < obs_min_dist){
           min_obs_index = nearby_lms[l].id;
           obs_min_dist = obs_to_lm;
         }
       }

       // Getting the nearest landmarks position
       for (unsigned int m = 0; m < nearby_lms.size(); m++){
         if (nearby_lms[m].id == min_obs_index){
           mu_x = nearby_lms[m].x;
           mu_y = nearby_lms[m].y;
           break;
         }
       }

       expon = pow(x_m-mu_x,2)/(2*pow(sig_x,2))+ pow(y_m-mu_y,2)/(2*pow(sig_y,2));
       weight *= gauss_norm * exp(-expon);
      }
      particles[i].weight = weight;
      weights[i] = weight;
   }	
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   vector<Particle> resampled_particles (num_particles);

   double beta = 0;
   int ind = rand() % num_particles;
   double w_max = *max_element(weights.begin(), weights.end());

   for (int i = 0; i < num_particles; i++){
     beta = beta * (rand() / (RAND_MAX + 1.0)) * (2*w_max);
     while(weights[ind]<beta){
            beta -= weights[ind];
            ind = (ind+1) % num_particles;
     }
     resampled_particles[i] = particles[ind];
   }
   particles = resampled_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
