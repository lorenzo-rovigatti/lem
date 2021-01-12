#!/usr/bin/env python3

import baggianalysis as ba
import lem

top_parser = ba.oxDNA_topology.TSP("topology.dat")
parser = ba.OxDNAParser(top_parser)
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file("trajectory.dat")

#ch = lem.ConvexHullAnalysis(trajectory)
#ch.print_to_file()

arms = top_parser.N_arms(0)
N_per_arm = top_parser.N_monomers_per_arm(0)
# this is the central bead
id_centre = 0
# these are the tips of the arms
id_particles = [a * N_per_arm + N_per_arm for a in range(arms)]
la = lem.LocalAnalysis(trajectory, id_centre=id_centre, id_particles=id_particles, find_correspondence=True, relative_to_centre=True)
la.print_to_file(prefix="lem_centre_")
# la = lem.LocalAnalysis(trajectory, id_centre=id_centre, id_particles=id_particles, find_correspondence=True)
# la.print_to_file(prefix="lem_centre_")
