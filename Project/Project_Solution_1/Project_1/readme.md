# Homework5-Final Project

## 1. file structure

All the task relevant head file are located in src/Topp/include/

1. cubic_spline_banded_test.h  -> smooth path generation using banded system
2. polymap.hpp -> generate the polyhedron obstacles and provide the signed distance and gradient
3. topp.h  ->  Time optimal trajectory optimization

  Program entrance: cubic_test.cpp  located in src/Topp/src/

noted: some support head file are located in src/Topp/include/support/



## 2. usage

```
cd Joeyyu-Project5
catkin build (catkin_make also works)
source devel/setup.zsh (or source devel/setup.bash)
roslaunch topp cubic_test.launch 
```

it will open the rviz and two rqt-plot window.

utilizing 2d Nav Goal in rviz to select the starts and goals (there is no visualization of the start and goals, you can choose any point around the obstacles)

wait for several seconds, it will display the smooth blue color path first.

After the topp optimization finished,  the rqt-plot will start plotting the velocity and acceleration. And the rviz will display the movement of the trajectory represented by a red ball.





