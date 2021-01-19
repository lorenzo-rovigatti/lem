#!/usr/bin/env python3
'''
Here we analyse LAMMPS trajectories of a single star polymer, taking the tips of the arms
as points for the analysis
'''
import baggianalysis as ba
import lempy
import sys

if len(sys.argv) < 4:
    print("Usage is %s trajectory N_arms N_per_arm" % sys.argv[0], file=sys.stderr)
    exit(1)
    
trajectory_file = sys.argv[1]
arms = int(sys.argv[2])
N_per_arm = int(sys.argv[3])

id_particles = [(a + 1) * N_per_arm + 1 for a in range(arms)]
id_centre = 1

parser = ba.LAMMPSDumpParser()
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file(trajectory_file)

ch = lempy.ConvexHullAnalysis(trajectory)
ch.print_to_file()

la = lempy.LocalAnalysis(trajectory, id_centre=id_centre, id_particles=id_particles, find_correspondence=True, relative_to_centre=True)
la.print_to_file()
