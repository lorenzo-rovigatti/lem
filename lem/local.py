'''
Created on 11 gen 2021

@author: lorenzo
'''

import numpy as np
from scipy.spatial import distance_matrix
from scipy.optimize import linear_sum_assignment

from lem.icp import icp
from lem.utils import make_pmf


class LocalAnalysis():

    def __init__(self, trajectory, id_centre=None, id_particles=None, do_linear_assigment=True):
        self.trajectory = trajectory
        self.id_centre = id_centre
        self.id_particles = id_particles
        self.do_linear_assignment = do_linear_assigment
        
        if id_particles is None:
            self.trajectory.reset()
            system = self.trajectory.next_frame()
            self.N = system.N()
        else:
            self.N = len(id_particles)
            
        self._compute_reference_matrix()
        self._analyse()
        
    def _extract_positions(self, system):
        centre = np.array([0., 0., 0.])
        if self.id_particles is None and self.id_centre is None:
            positions = system.positions()
        else:
            positions = []
            for p in system.particles():
                if self.id_centre is not None and p.index == self.id_centre:
                    centre = p.position
                elif self.id_particles is None or p.index in self.id_particles:
                    positions.append(p.position)
        
        return positions - centre
    
    def _assign_permutation(self, positions, base_positions):
        if self.do_linear_assignment:
            # compute the distance matrix
            graph_weights = distance_matrix(base_positions, positions)
            # perform the linear assignment
            row_ind, col_ind = linear_sum_assignment(graph_weights)
            # reorder the array according to the result of the linear assignment
            return positions[col_ind,:]
        else:
            return positions

    def _compute_reference_matrix(self):
        self.trajectory.reset()
        system = self.trajectory.next_frame()
        n_confs = 0
        base_positions = None
        while system != None:
            n_confs += 1
            
            positions = self._extract_positions(system)
            
            if base_positions is None:
                base_positions = np.copy(positions)
                self.reference_positions = np.copy(base_positions)
            else:
                T, distances, iterations = icp(positions, base_positions)
        
                # make C a homogeneous representation of src
                C = np.ones((self.N, 4))
                C[:, 0:3] = np.copy(positions)
                # transform C according to the change of coordinates obtained with the icp procedure
                C = np.dot(T, C.T).T
                # obtain the new coordinates
                positions_corrected = C[:,:3]
                
                self.reference_positions += self._assign_permutation(positions_corrected, base_positions)
            
            system = self.trajectory.next_frame()
        
        self.reference_positions /= n_confs
        
        dest_distance_matrix = distance_matrix(self.reference_positions, self.reference_positions)
        cut_off = dest_distance_matrix.max()
        dest_distance_matrix_norm = dest_distance_matrix / cut_off
        
        first_piece = dest_distance_matrix_norm < 0.5
        second_piece = (dest_distance_matrix_norm >= 0.5) & (dest_distance_matrix_norm < 1.0)
        third_piece = dest_distance_matrix_norm >= 1.0
        
        weight_matrix = np.piecewise(
            dest_distance_matrix_norm,
            [first_piece, second_piece, third_piece],
            [lambda r: 1 - 6 * r ** 2 + 6 * r ** 3, lambda r: 2 - 6 * r + 6 * r ** 2 - 2 * r ** 3, 0.]
        )
        # see https://stackoverflow.com/a/54637261/5140209
        self.weight_matrix = np.broadcast_to(weight_matrix[:,:, np.newaxis, np.newaxis], (self.N, self.N, 3, 3))
        
        # now we compute D_inv for each reference point
        self.reference_diff_matrix = self.reference_positions[:, np.newaxis,:] - self.reference_positions
        # first we perform an outer product on the last axis
        D = self.reference_diff_matrix[:,:,:, np.newaxis] * self.reference_diff_matrix[:,:, np.newaxis,:] * self.weight_matrix
        # then we sum over all the points (i.e. on the second axis)
        D = np.sum(D, axis=1)
        # and invert the matrix
        self.D_inv = np.linalg.inv(D)
        
    def _compute_F(self, system):
        positions = self._extract_positions(system)
        
        T, distances, iterations = icp(positions, self.reference_positions)
    
        # make C a homogeneous representation of src
        C = np.ones((self.N, 4))
        C[:, 0:3] = np.copy(positions)
        # transform C according to the change of coordinates obtained with the icp procedure
        C = np.dot(T, C.T).T
        # obtain the new coordinates
        positions_corrected = C[:,:3]
        
        positions_final = self._assign_permutation(positions_corrected, self.reference_positions)
        
        # compute A for each point
        diff_matrix = positions_final[:, np.newaxis,:] - positions_final
        A = diff_matrix[:,:,:, np.newaxis] * self.reference_diff_matrix[:,:, np.newaxis,:] * self.weight_matrix
        # then we sum over all the points (i.e. along the second axis)
        A = np.sum(A, axis=1)
    
        # return F
        return A @ self.D_inv
    
    def _analyse(self):
        self.J_all = []
        self.I_all = []
        self.J_conf = []
        self.I_conf = []
        
        self.trajectory.reset()
        system = self.trajectory.next_frame()
        confs = 0
        while system != None:
            confs += 1
            if confs % 100 == 0:
                print("Analysed %s confs" % confs)
        
            # and now we can compute the Cauchy-Green strain tensor
            F = self._compute_F(system)
            F_T = np.transpose(F, axes=(0, 2, 1))
            C = np.matmul(F_T, F)
        
            # from the invariants of C we obtain J and I
            J = np.sqrt(np.linalg.det(C))
            I = np.trace(C, axis1=1, axis2=2) / J ** (2. / 3.)
            
            self.J_all.append(J)
            self.I_all.append(I)
                
            F_avg = np.average(F, axis=0)
            C_avg = F_avg.T @ F_avg
            J_from_F_avg = np.sqrt(np.linalg.det(C_avg))
            self.J_conf.append(J_from_F_avg)
            self.I_conf.append(np.trace(C_avg) / J_from_F_avg ** (2. / 3.))
        
            system = self.trajectory.next_frame()
            
    def print_to_file(self, prefix="lem_"):
        with open("%sJ_all.dat" % prefix, "w") as J_out, open("%sI_all.dat" % prefix, "w") as I_out:
            for i in range(len(self.J_all)):
                np.savetxt(J_out, self.J_all[i])
                np.savetxt(I_out, self.I_all[i])
                # add an empty line after each set of results
                print("", file=J_out)
                print("", file=I_out)
            
        np.savetxt("%sJ_conf.dat" % prefix, self.J_conf)
        np.savetxt("%sI_conf.dat" % prefix, self.I_conf)

        pmf_J_conf = make_pmf(self.J_conf)
        pmf_I_conf = make_pmf(self.I_conf)

        np.savetxt("%sJ_conf_pmf.dat" % prefix, pmf_J_conf)
        np.savetxt("%sI_conf_pmf.dat" % prefix, pmf_I_conf)
    
