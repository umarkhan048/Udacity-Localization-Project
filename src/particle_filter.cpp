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

   for (unsigned int i=0; i < observations.size(); i++){
	double minDist;
	int mapID;
	for (unsigned int j=0; j < predicted.size(); j++){
		double distance = sqrt(pow((predicted[j].x - observations[i].x), 2) + pow((predicted[j].y - observations[i].y), 2));

	if (j==0){
		minDist = distance;
	}
	else if (distance < minDist) {
		minDist = distance;
		mapID = predicted[j].id;
	}
	}
	observations[i].id = mapID;
   }
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
   for (int i = 0; i < num_particles; i++){

	vector<LandmarkObs> par_inrange_lm;
	double x = particles[i].x;
	double y = particles[i].y;
	double theta = particles[i].theta;

	for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++){

	   double lmX = map_landmarks.landmark_list[j].x_f;
	   double lmY = map_landmarks.landmark_list[j].y_f;
	   int lmID = map_landmarks.landmark_list[j].id_i;
	   // Get landmark to particle distance;
	   double lm_dist = dist(x, lmX, y, lmY);

	   if (lm_dist <= sensor_range){
		par_inrange_lm.push_back(LandmarkObs{lmID, lmX, lmY});
	   }
	}

	// Transform from vehicle to particle coordinates.
    	vector<LandmarkObs> transformed_obs;
    	for(unsigned int j = 0; j < observations.size(); j++) {
           double transformed_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
           double transformed_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
           transformed_obs.push_back(LandmarkObs{ observations[j].id, transformed_x, transformed_y });
    	}

	// Associating observations to landmarks by nearest neighbor method
	dataAssociation(par_inrange_lm, transformed_obs);

	// Updating weights
	for (unsigned int k = 0; k < transformed_obs.size(); k++){
	   double obs_x = transformed_obs[k].x;
	   double obs_y = transformed_obs[k].y;
	   double mu_x = par_inrange_lm[k].x;
	   double mu_y = par_inrange_lm[k].y;

	   double gauss_norm;
	   gauss_norm = 1/(2 * M_PI * std_landmark[0] * std_landmark[1]);

	   double expon;
	   expon = (pow(obs_x - mu_x,2) / (2*pow(std_landmark[0],2))) + (pow(obs_y - mu_y,2) / (2*pow(std_landmark[1],2)));

	   double weight;
	   weight = gauss_norm * exp(expon);
	   particles[i].weight = weight;
	}

   }	
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   // Get weights and max weight.

   double max_weight;
   vector<double> par_weights;

   // find max weight
   for (unsigned int i = 0; i < num_particles; i++){
	par_weights.push_back(particles[i].weight);

	if (particles[i].weight > max_weight){
	   max_weight = particles[i].weight;
	}
   }

   double beta = 0.0;
   std::default_random_engine gen;
   // for 2 * random(0-> max_weight) * max_weight 
   std::uniform_real_distribution<double> random_double(0.0, max_weight);

   // generate indexes
   std::uniform_int_distribution<int> index_distrib(0, num_particles-1);

   int index = index_distrib(gen);

   // vector for resampled particles
   vector<Particle> resamp_particles;
 
   for (unsigned int j = 0; j < num_particles; j++){
	beta = beta + random_double(gen) * 2.0;
	while (beta > weights[index]){
	   beta -= par_weights[index];
	   index = (index + 1) % num_particles;
	}
	resamp_particles.push_back(particles[index]);
   }
   particles = resamp_particles; 
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
