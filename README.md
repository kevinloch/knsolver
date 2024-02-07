# knsolver
knsolver is a tool for exploring/analyzing the Kerr-Newman metric near Planck scale

It runs in the [Julia Language](https://julialang.org/) interpreter and renders plots with [gnuplot](http://www.gnuplot.info/)

## Usage
- Edit knsolver.jl to set/change configuration parameters
- Start Julia interpreter
- Change to working directory

  `julia> cd("~/knsolver/src")`
- Load knsolver

  `julia> include("knsolver.jl")`
- run

  `julia> main()`

- knsolver will output a CSV file 'tmp.csv' and attempt to run gnuplot to generate 'tmp.pdf'. On Mac OSX it will also attempt to load 'tmp.pdf' in the Preview app. You can edit this behavior in knsolver.jl to suit your workflow.
- If you want to make parameter changes, edit knsolver.jl again, and repeat the include and main() steps.

## Background
The inspiration for this tool was to explore *neutral metric* solutions to General Relativity (where dτ/dt=1), but it can be used to plot solutions to any target time ratio, including zero. That non-trivial, finite, real, neutral metric solutions exist may be surprising even to physicists, who are used to focusing on effects occuring near the zero metric. Indeed, in the trivial case with M != 0, Q = 0, J = 0, and zero test particle velocity there are no finite solutions for dτ/dt=1. However, if M and any one of Q, J, or test particle velocity is non-zero, one or more solutions of finite radius exist.

## Negative Mass?
General relativity is time-symmetric. The KN metric for example is specified in terms of time squared, so that *identical* positive (foward time) and negative (reverse time) solutions exist for any given set of parameters. It is *not* however symmetric in +/- Mass. This due to having mass terms (r_s and a) that contain M that are not quadratic. If we disregard negative mass (energy) as being unphysical we are still left with positive mass operating in reverse time, which behaves *exactly* like negative mass (both gravitational and inertial) would in forward time. By combining these two features (symmetric time behavior and asymmetric mass behavior) we can accept negative values for M, understanding that it just represents positive mass operating in reverse time.

One interesting consequence is that while some solutions to the neutral metric give negative values for r (which we also consider unphysical), inverting the sign of M also inverts the sign of r, turning a negative radius into a positive radius.

## This may lead to an understanding of how Standard Model (SM) particles gain mass
(or the value of the Higgs field couplings in physicist speak)
As can be shown by knsolver plots, very nonlinear and interesting things happen right around Planck scale. If a pair of positive and negative Planck scale masses were able to coexist in a stable resonance (with some gravitational "partitioning" and perhaps enforced by time neutrality) and the positive mass was slightly larger than the negative mass, there would be a tiny residual positive mass in the far field. Some version of this model might explain the [Hierarchy problem](https://en.wikipedia.org/wiki/Hierarchy_problem), (why SM particle masses so tiny compared to the Planck mass).

### What about antiparticles?
There is nothing in this model that requires antiparticles to have negative mass, or to propogate in negative time. In fact antiparticles are known to propogate in forward time, and are known to have positive mass. The preference for forward time in macroscopic (>> Planck scale) behavior could be an intrinsic result of the non-linear solutions regardless of the chosen parameters.
