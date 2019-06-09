# RandomDynamics

A toy projects that shows how simple rules can lead to interesting results and complex behaviour. The idea is inspired by cellular automata. In contrast to cellular automata, this approach does not use a grid. In order to calculate the interaction of individual agents more efficiently, a tree algorithm was used.

The rule how these agents behave is really simple. Every agent just adapts the behaviour (direction of movement) from other agents that are within a distance `DIST`.

Starting with a random initial state of agents that looks like the following figure

<div align="center">
<img src="https://github.com/KaiFabi/RandomDynamics/blob/master/init.png" height="320">
</div>

interesting patterns like vortices emerge.

Here are some results for `DIST=0.03`, `DIST=0.07` and `DIST=0.14` (from left to right).

<p align="center">
![](https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p14.gif)
![](https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p07.gif)
![](https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p03.gif)
</p>
