## Unscented Kalman Filter

### Scope

The goal of this project is to track a bicycle or other target using a more accurate filter than the Extended Kalman Filter.  This Unscented Kalman Filter fuses lidar and radar sensor data to track a 5-dimensional state consisting of x and y position, speed (velocity magnitude), yaw angle, and yaw angle rate of change.  The model used is the _Constant Turn Rate and Velocity Magnitude Model_ (CTRV).

### Process

For both lidar and radar, the state and covariance are predicted using the UKF process.  Longitudinal and yaw acceleration are added as noise in the prediction.  Following prediction, a measurement update is applied.  Predicted sigma points are converted into the appropriate measurement space (lidar or radar) and the predicted measurement mean and covariance are computed.  Using the difference between measurement and prediction, the filter's state and covariance belief are then updated.

### Unscented vs Extended

Whereas the Extended Kalman Filter handles non-linear process models by linearization, the Unscented Kalman Filter generates a set of _sigma points_ representing a Gaussian distribution.  These points are then passed through the non-linear process function to create a predicted state distribution.  The resulting Gaussian distribution more accurately represents the actual, non-Gaussian, distribution.

### Accuracy

Following prediction and measurement update, the filter's state is compared to the _ground truth_ position and velocity.  The 4-dimensional position and velocity vector is computed from the 5-dimensional filter state and used to generate the RMSE.

### Consistency
A major component of this project involves tuning the process noise standard deviation parameters.  In addition to overall accuracy, _consistency_ is also evaluated for laser and radar samples.  The _Normalized Innovation Squared_ (NIS) values are plotted below.  These values should have a Chi-Squared distribution.  Statistically, for laser (2 degrees of freedom) 5% of values should be higher than 5.991.  For radar (3 degrees of freedom) 5% of values should be higher than 7.815.  The graphs below show that computed NIS values are visually consistent with the expected distributions.  Note that the first value for radar (69.64) is omitted to avoid scaling the graph based on initial noise.

If the NIS values are not distributed as expected, the process noise standard deviations are either overly conservative (too high; low NIS) or overly optimistic (too low; high NIS).

#### Laser
![alt text][nis_laser]
[nis_laser]: images/nis_laser.jpg "NIS Laser"

#### Radar
![alt text][nis_radar]
[nis_radar]: images/nis_radar.jpg "NIS Radar"

### Difficulties

This project required tuning the initial covariance matrix and the acceleration standard deviations, and the RMSE could therefore be expected to start out poor.  In addition to this, angle normalization must be applied when a difference of angles is computed.  I made the mistake of questioning this requirement for trigonometry reasons, and only realized later that this logic probably does not apply to covariance computations.  When I first went down this path, I normalized angles in multiple places, but limiting to points at which a difference was computed produced better results.  During all of this I was fighting numeric instability and a probable failure to initialize a variable somewhere.  It finally turned out I had missed one MatrixXd initialization.  Before I spotted that bug, I had gotten close to the required RMSE but couldn't quite hit the mark.  Parameters were very touchy and small changes often caused dramatic shifts in RMSE.


