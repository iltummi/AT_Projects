# Coupled GPS/INS integration for UAV Applications
Sensor fusion algorithm with extended kalman filter; the first EKF is used to estimate drone's attitude and gyroscope bias thanks to attitude measure obtained with Q-TRIAD algorithm, while the other EKF is used to estimate drone's positon, velocity and acceloremeter bias. GPS and IMU sensors are simlauted thanks to MATLAB's gpsSensor and imuSensor function, avaiable in the Navigation Toolbox.
Desidered trajectory is a circle around a fixed coordinate and during this path I supposed a sinusoidal attitude with different amplitude along yaw, pitch and roll; this trajectory is simulated with waypointTrajectory function, in the same toolbox.
To run the algorithm just run the main file.
