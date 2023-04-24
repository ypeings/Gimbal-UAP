# Gimbal-UAP
A Python code to examine potential flight paths for the Gimbal UAP

1/ The user needs to set a few parameters for the encounter at the top of the code: 
- wind direction (angle in degrees, relative to initial heading) and speed in Knots ("the wind is 120 Knots to the West")
- the type of flight path for the object : const_head or const_alt. 
- const_head plots a straight line through the lines of sight (LOS), starting at an initial distance, in Nautical miles Nm, set by the "dist_init" parameter
  The direction of flight is specified by the "offset" parameter, as an angle of deviation from wind direction in degrees ("they're going against the wind")
- const_alt plots a leveled flight path at a certain altitude through the LOS ("Alt_object" parameter, altitude in feet)
- the field of view (FOV) is set to 0.35° (NAR2 ATFLIR FOV)
- the smth_LOS parameter allows to choose for the option of smoothing the LOS using angular motion of the clouds (True by default)

2/ The code plots the F-18 flight path for the encounter, using indicated air speed (converted to true air speed) and bank angle (converted to rate of turn).
   The F-18 flight path (with or without the effect of selected wind) is plotted in Figures/Fig_F18_path.pdf

3/ The code then plots the associated LOS, and the UAP flight path that corresponds to the set of parameters specified by the user. 
   An interactive window appears that shows both flight paths and the LOS, in a 3D space (X, Y, Z coordinates in Nm). The persepctive, zoom, can be changed 
   to visualize the results, and plots for the encounter can be saved. 
   Several graphs are generated that indicate air speed and ground speed for the object (Fig_speed.pdf), altitude, (Fig_alt.pdf) and other parameters, 
   in function of time (frame #) 


The default flight path obeys the following conditions that align with the context of the encounter:
- 6-8Nm range -> starts at 8Nm
- 15% increase in size throughout the video-> -15% decrease in distance
- straight flight path, stopping/reversing direction (going throught the LOS)
- object roughly facing the wind (towards the East)
- object going left to right on the SA
- stern conversion (F-18 originally going behind the object)
