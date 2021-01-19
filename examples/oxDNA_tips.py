#!/usr/bin/env python3
'''
Here we analyse oxDNA trajectories of a single star polymer, taking the tips of the arms
as points for the analysis
'''
import baggianalysis as ba
import lempy
import sys

if len(sys.argv) < 3:
    print("Usage is %s topology trajectory" % sys.argv[0], file=sys.stderr)
    exit(1)
    
topology_file = sys.argv[1]
trajectory_file = sys.argv[2]

top_parser = ba.oxDNA_topology.TSP(topology_file)
parser = ba.OxDNAParser(top_parser)
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file(trajectory_file)

# analyse with the convex hull
ch = lempy.ConvexHullAnalysis(trajectory)
ch.print_to_file()

arms = top_parser.N_arms(0)
N_per_arm = top_parser.N_monomers_per_arm(0)
# in an oxDNA configuration the central bead has always index 0
id_centre = 0
# these are the tips of the arms
id_particles = [a * N_per_arm + N_per_arm for a in range(arms)]
la = lempy.LocalAnalysis(trajectory, id_centre=id_centre, id_particles=id_particles, find_correspondence=True, relative_to_centre=True)
la.print_to_file()
