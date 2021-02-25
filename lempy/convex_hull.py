import numpy as np
import baggianalysis as ba

from .utils import make_pmf


class ConvexHullAnalysis():

    def __init__(self, trajectory, also_spherical=False):
        self.also_spherical = also_spherical
        self._analyse(trajectory)
        
    def _analyse(self, trajectory):
        ch = ba.ConvexHull()
        
        trajectory.reset()
        system = trajectory.next_frame()
        self.eigenvalues = []
        while system != None:
            ch.analyse_system(system)
            self.eigenvalues.append(self._eigenvalues_from_ch_results(ch.result()))
            system = trajectory.next_frame()

        self.V = 4. * np.pi * np.sqrt(3) * np.prod(np.sqrt(self.eigenvalues), axis=1)
        self.V_avg = np.average(self.V)

        R2 = np.average(self.eigenvalues, axis=0)
        strains = self.eigenvalues / R2

        self.J = np.sqrt(np.prod(strains, axis=1))
        self.I = np.sum(strains, axis=1) * self.J ** (-2. / 3.)
        
        if self.also_spherical:
            R2 = np.average(self.eigenvalues)
            strains = self.eigenvalues / R2
    
            self.J_spherical = np.sqrt(np.prod(strains, axis=1))
            self.I_spherical = np.sum(strains, axis=1) * self.J ** (-2. / 3.)

    def print_to_file(self, prefix="ch_"):
        np.savetxt("%sJ.dat" % prefix, self.J)
        np.savetxt("%sI.dat" % prefix, self.I)
        np.savetxt("%sV.dat" % prefix, self.V)

        pmf_J = make_pmf(self.J)
        pmf_I = make_pmf(self.I)

        np.savetxt("%sJ_pmf.dat" % prefix, pmf_J)
        np.savetxt("%sI_pmf.dat" % prefix, pmf_I)
        
        if self.also_spherical:
            np.savetxt("%sJ_spherical.dat" % prefix, self.J_spherical)
            np.savetxt("%sI_spherical.dat" % prefix, self.I_spherical)
            
            pmf_J = make_pmf(self.J_spherical)
            pmf_I = make_pmf(self.I_spherical)
    
            np.savetxt("%sJ_spherical_pmf.dat" % prefix, pmf_J)
            np.savetxt("%sI_spherical_pmf.dat" % prefix, pmf_I)

    def _eigenvalues_from_ch_results(self, ch_result):
        triangle_coms = np.array(list(map(lambda t: (t.v1 + t.v2 + t.v3) / 3., ch_result.triangles)))
        ch_com = np.average(triangle_coms, axis=0)
        triangle_coms -= ch_com
        
        # https://stackoverflow.com/questions/62153830/how-do-i-efficiently-compute-the-gyration-tensor-in-numpy
        gyr_tensor = np.einsum('im,in->mn', triangle_coms, triangle_coms) / triangle_coms.shape[0]
        ev, _ = np.linalg.eig(gyr_tensor)
        ev = sorted(ev)
    
        return ev
