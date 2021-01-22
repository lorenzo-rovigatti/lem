'''
Created on 11 gen 2021

@author: lorenzo
'''

import numpy as np
from scipy.spatial import distance_matrix
from scipy.optimize import linear_sum_assignment

from .icp import icp, best_fit_transform
from .utils import make_pmf


class LocalAnalysis():

    def __init__(self, trajectory, id_centre=None, id_particles=None, find_correspondence=True, relative_to_centre=False):
        self.trajectory = trajectory
        self.id_centre = id_centre
        self.id_particles = id_particles
        self.find_correspondence = find_correspondence
        if self.find_correspondence:
            self.align_function = icp
        else:
            self.align_function = best_fit_transform
        self.relative_to_centre = relative_to_centre
        
        if id_particles is None:
            self.trajectory.reset()
            system = self.trajectory.next_frame()
            self.N = system.N()
            if id_centre is not None:
                self.N -= 1
            self.trajectory.reset()
        else:
            self.N = len(id_particles)
            
        self._compute_reference_matrix()
        self._analyse()
        
    def _extract_positions(self, system):
        if self.id_particles is None and self.id_centre is None:
            positions = system.positions()
        else:
            positions = []
            for p in system.particles():
                if self.id_centre is not None and p.index == self.id_centre:
                    centre = p.position
                elif self.id_particles is None or p.index in self.id_particles:
                    positions.append(p.position)
            
        # we compute the average ourselves        
        if self.id_centre is None:
            centre = np.average(positions, axis=0)
        
        return positions - centre
    
    def _align(self, to_be_aligned, reference):
        T = self.align_function(to_be_aligned, reference)[0]
        
        if self.relative_to_centre:
            rot_matrix = T[:3,:3]
            C = np.copy(to_be_aligned)
            return np.dot(rot_matrix, C.T).T
        else:
            # make C a homogeneous representation of to_be_aligned
            C = np.ones((self.N, 4))
            C[:, 0:3] = np.copy(to_be_aligned)
            # transform C according to the change of coordinates obtained with the alignment procedure
            C = np.dot(T, C.T).T
            # return the new coordinates
            return C[:,:3]
    
    def _assign_permutation(self, positions, base_positions):
        if self.find_correspondence:
            # compute the distance matrix
            graph_weights = distance_matrix(base_positions, positions)
            # perform the linear assignment
            _, col_ind = linear_sum_assignment(graph_weights)
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
                positions_corrected = self._align(positions, base_positions)
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
        
        ref_pos_squared = self.reference_positions * self.reference_positions
        self.D_global = np.broadcast_to(ref_pos_squared[:,np.newaxis,:], (self.N, 3, 3))
        self.D_global = np.sum(self.D_global, axis=0)
        
    def _compute_F(self, system):
        positions = self._extract_positions(system)
        positions_corrected = self._align(positions, self.reference_positions)
        positions_final = self._assign_permutation(positions_corrected, self.reference_positions)
        
        if self.relative_to_centre:
            A = np.sum(positions_final[:,:,np.newaxis] * self.reference_positions[:,np.newaxis,:], axis=0)
            
            F_global = A / self.D_global
        else:
            F_global = None
            
        # compute A for each point
        diff_matrix = positions_final[:, np.newaxis,:] - positions_final
        A = diff_matrix[:,:,:, np.newaxis] * self.reference_diff_matrix[:,:, np.newaxis,:] * self.weight_matrix
        # then we sum over all the points (i.e. along the second axis)
        A = np.sum(A, axis=1)
    
        # F = A * D^-1
        return A @ self.D_inv, F_global
    
    def _analyse(self):
        self.J_local = []
        self.I_local = []
        self.J_local_avg = []
        self.I_local_avg = []
        self.J_global = []
        self.I_global = []
        
        self.trajectory.reset()
        system = self.trajectory.next_frame()
        confs = 0
        while system != None:
            confs += 1
            if confs % 100 == 0:
                print("Analysed %s confs" % confs)
        
            # and now we can compute the Cauchy-Green strain tensor
            F, F_global = self._compute_F(system)
            F_T = np.transpose(F, axes=(0, 2, 1))
            C = np.matmul(F_T, F)
        
            # from the invariants of C we obtain J and I
            J = np.sqrt(np.linalg.det(C))
            I = np.trace(C, axis1=1, axis2=2) / J ** (2. / 3.)
            
            self.J_local.append(J)
            self.I_local.append(I)
            
            F_avg = np.average(F, axis=0)
            C_avg = F_avg.T @ F_avg
            J_from_F_avg = np.sqrt(np.linalg.det(C_avg))
            self.J_local_avg.append(J_from_F_avg)
            self.I_local_avg.append(np.trace(C_avg) / J_from_F_avg ** (2. / 3.))

            if F_global is not None:                
                C_global = F_global.T @ F_global
                J_from_F_global = np.sqrt(np.linalg.det(C_global))
                self.J_global.append(J_from_F_global)
                self.I_global.append(np.trace(C_global) / J_from_F_global ** (2. / 3.))
        
            system = self.trajectory.next_frame()
            
    def print_to_file(self, prefix="lem_"):
        with open("%sJ_local.dat" % prefix, "w") as J_out, open("%sI_local.dat" % prefix, "w") as I_out:
            for i in range(len(self.J_local)):
                np.savetxt(J_out, self.J_local[i])
                np.savetxt(I_out, self.I_local[i])
                # add an empty line after each set of results
                print("", file=J_out)
                print("", file=I_out)
                
        np.savetxt("%sJ_local_avg.dat" % prefix, self.J_local_avg)
        np.savetxt("%sI_local_avg.dat" % prefix, self.I_local_avg)
                
        pmf_J_local_avg = make_pmf(self.J_local_avg)
        pmf_I_local_avg = make_pmf(self.I_local_avg)
        
        np.savetxt("%sJ_local_avg_pmf.dat" % prefix, pmf_J_local_avg)
        np.savetxt("%sI_local_avg_pmf.dat" % prefix, pmf_I_local_avg)
            
        if self.relative_to_centre:
            np.savetxt("%sJ_global.dat" % prefix, self.J_global)
            np.savetxt("%sI_global.dat" % prefix, self.I_global)
    
            pmf_J_global = make_pmf(self.J_global)
            pmf_I_global = make_pmf(self.I_global)
    
            np.savetxt("%sJ_global_pmf.dat" % prefix, pmf_J_global)
            np.savetxt("%sI_global_pmf.dat" % prefix, pmf_I_global)
    
