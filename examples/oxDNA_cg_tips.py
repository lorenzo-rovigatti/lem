#!/usr/bin/env python3
'''
Here we analyse oxDNA trajectories of a single star polymer, taking the average position of the 
last X monomers that compose each arm as points for the analysis
'''
import baggianalysis as ba
import lempy
import sys

if len(sys.argv) < 4:
    print("Usage is %s topology trajectory number_of_monomers_to_average_over" % sys.argv[0], file=sys.stderr)
    exit(1)
    
topology_file = sys.argv[1]
trajectory_file = sys.argv[2]
average_over = int(sys.argv[3])
if average_over < 2:
    print("number_of_monomers_to_average_over should be an integer number larger than 1", file=sys.stderr)
    exit(1)

top_parser = ba.oxDNA_topology.TSP(topology_file)
# in an oxDNA configuration the central bead has always index 0
id_centre = 0
arms = top_parser.N_arms(0)
N_per_arm = top_parser.N_monomers_per_arm(0)

# we now build lists of ids over which we average to get each arm's "coarse-grained" position
# the first list contains only the central bead
id_lists = [[id_centre,],]
for a in range(arms):
    new_list = [(a + 1) * N_per_arm - monomer for monomer in range(average_over)]
    id_lists.append(new_list)
    
my_filter = ba.MapParticles(id_lists)
parser = ba.OxDNAParser(top_parser)
trajectory = ba.FullTrajectory(parser)
trajectory.add_filter(my_filter)
trajectory.initialise_from_trajectory_file(trajectory_file)

# analyse with the convex hull
ch = lempy.ConvexHullAnalysis(trajectory)
ch.print_to_file("ch_cg%d_" % average_over)

la = lempy.LocalAnalysis(trajectory, id_centre=id_centre, find_correspondence=True, relative_to_centre=True)
la.print_to_file("lem_cg%d_" % average_over)
