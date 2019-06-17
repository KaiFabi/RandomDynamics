# RandomDynamics

This toy project shows how simple rules can lead to interesting results and complex behaviour. The idea is inspired by cellular automata. In contrast to cellular automata, this approach does not use a grid. In order to calculate the interaction of individual agents more efficiently, a tree algorithm was used.

The rule that determines how these agents behave is quite simple. Each agent adapts his movement to that of the surrounding agents within distance `DIST`.

Starting with a random initial state of agents that looks like the following figure

<div align="center">
<img src="https://github.com/KaiFabi/RandomDynamics/blob/master/init.png" height="500">
</div>

interesting patterns like vortices emerge.

Here are some results for `DIST=0.03`, `DIST=0.07` and `DIST=0.14` (from left to right).

<p align="center">
<img src="https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p14.gif" height="300">
<img src="https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p07.gif" height="300">
<img src="https://github.com/KaiFabi/RandomDynamics/blob/master/output_dist_0p03.gif" height="300">
</p>

Compile and run the program using

`gcc -O -Wall automata.c -o automata -lm`
and 
`./automata`

$`\sqrt{2}`$
