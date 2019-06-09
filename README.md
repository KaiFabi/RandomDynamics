# RandomDynamics

A toy projects that shows how simple rules can lead to interesting results and complex behaviour. The idea is inspired by cellular automata. In contrast to cellular automata, this approach does not use a grid. In order to calculate the interaction of individual agents more efficiently, a tree algorithm was used.

The rule how these agents behave is really simple. Every agent just adapts the behaviour (direction of movement) from other agents that are within a distance `DIST`.

Here are some results for different distances:

![](https://github.com/KaiFabi/RandomDynamics/blob/master/ouput_dist_0p14.gif)
![](https://github.com/KaiFabi/RandomDynamics/blob/master/ouput_dist_0p07.gif)
![](https://github.com/KaiFabi/RandomDynamics/blob/master/ouput_dist_0p03.gif)
