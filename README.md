# Trajectory planning for a spacecraft from Earth to Enceladus

The purpose of the project is to simulate an interplanetary mission from Earth to 
Enceladus, a satellite of Saturn. This was accomplished by using the approximation 
given by the "Patched Conics" method, which allowed the mission to be divided into a series of phases, in each of which the spacecraft is subjected primarily under 
the gravitational action of a specific celestial body. The resolution of the 
two-body problem, applied to each phase of the mission, allows us to uniquely obtain 
a closed form of the trajectory in space of the secondary body (spacecraft), 
orbiting the main body (celestial body). 
The conic arcs of the trajectory associated with each phase will need to be
appropriately connected. This requires the use of additional simplifying assumptions; 
it is necessary to introduce the concept of a sphere of influence.
