import baggianalysis as ba
import numpy as np
import lem

parser = ba.GenericOxDNAParser("ba_topology.dat")
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file("trajectory.dat")

ch = lem.ConvexHullAnalysis(trajectory)
ch.print_to_file()

