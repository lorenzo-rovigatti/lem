import numpy as np

def make_pmf(x, remove_last=None):
    hist, bins = np.histogram(x, bins='auto', density=True)
    bins = bins[:-1] + (bins[1] - bins[0]) / 2.
    pmf = np.column_stack((bins, hist))
    # get rid of the last points
    if remove_last is not None:
        if len(bins) < 2 * remove_last:
            remove_last = len(bins) / 5
        pmf = pmf[pmf[:,1] > 0][:-remove_last]
    else:
        pmf = pmf[pmf[:,1] > 0]
    pmf[:,1] = -np.log(pmf[:,1])

    return pmf

def print_data_and_pmf_to_file(data, prefix):
    np.savetxt("%s.dat" % prefix, data)
    pmf = make_pmf(data)
    np.savetxt("%s_pmf.dat" % prefix, pmf)
    

def compute_elastic_moduli(J_par, I_par, V_avg):
    K = 2. * J_par / V_avg
    G = 2. * I_par / V_avg
    Y = 9. * K * G / (3. * K + G)
    nu = (3. * K - Y) / (6. * K)

    D = 1.5 * (1. - nu**2.) / Y
    hertzian_coeff = 1. / (5. * D)

    results = {
        'K' : K,
        'G' : G,
        'Y' : Y,
        'nu' : nu,
        'D' : D,
        'coeff' : hertzian_coeff
    }

    return results
