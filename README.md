# Overview
This repository contains the implementation of the localization project of Udacity's Self-Driving Car Nanodegree.

The flowchart below shows the steps of the particle filter algorithm as well as its inputs.

![Particle Filter Flowchart](/images/Particle_Filter_Algo_Flowchart.png)
Source: Udacity

### Initialization

In the initialization step, the total number of particles is set and the particles are initialized. The first position is set to the GPS position. Random Gaussian noise is added.

### Prediction

In the prediction step, the yaw rate and velocity information is used to predict the x and y and theta of each particle in the time period delta_t. The case of very low yaw rate is also considered. Gaussian noise is added to each particle's x, y and theta values to immitate sensor noise.

### Update Step (Update of Weights)

In the update step, first of all, the landmarks in the sensor range are selected. Then all the observations are transformed from the vehicle coordinates to particle coordinates system. In the next step, the association between nearby landmarks and observations is performed. For the purpose of this project, the data association is implemented based on the shortest distance between landmark and observation.

After above steps, all the information required to calculate particle's final weight is available. The particles final weight will be calculated as the product of each measurement's Multivariate-Gaussian probability density.

### Resample

In the resample step, the particles are resampled with probability proportioanl to their weight.
