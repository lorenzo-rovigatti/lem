# lem

This repo contains the `lempy` python package and a few accompanying scripts. `lempy` can analyse simulation trajectories to extract the local and global fluctuations of the shape and volume of soft objects such as star polymers or microgels.

The analysis can be carried out either with the convex hull construction or with a different method that is described in the `docs/theory.pdf` file and can be used to extract information about the *local* elastic fluctuations.

## Installation

The package (and its requirements) can be installed with 

```
git clone http://github.com/lorenzo-rovigatti/lem.git
cd lem
pip3 install .
```

You may need to use `pip` instead of `pip3` (or specify `pip`'s path), depending on your local machine's details.

## Usage

### Initialising a trajectory object

The basic data structure that should be fed to `lempy` is a simulation trajectory. This is done by using [`baggianalysis`](https://github.com/lorenzo-rovigatti/baggianalysis), which supports several trajectory formats. Here is an example that assumes a [LAMMPS](https://lammps.sandia.gov/) trajectory file:

```Python
import baggianalysis as ba
import lempy

parser = ba.LAMMPSDumpParser()
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file("trajectory.lammpstrj")
```

Once the trajectory has been initialised, `lempy` can be used to analyse it. 

### Convex hull analysis

The convex hull construction and analysis can be run with

```Python
ch = lempy.ConvexHullAnalysis(trajectory)
ch.print_to_file()
```

The `print_to_file` call will generate five files which will contain the instantaneous shape and volume fluctuations (`ch_I.dat` and `ch_J.dat`), their normalised histogram (`ch_I_pmf.dat` and `ch_J_pmf.dat`) and the istantaneous volume (`ch_V.dat`).

The constructor takes an optional parameter (`also_spherical`) which, if set to `True`, will print four additional files (`ch_I_spherical.dat`, `ch_J_spherical.dat`, `ch_I_spherical_pmf.dat` and `ch_J_spherical_pmf.dat`) containing the results obtained by using a sphere in place of an ellipsoid as the reference configuration. 

### LEM analysis

The local-elasticity analysis can be carried out in a similar fashion by using the `LocalAnalysis` object as follows:

```Python
la = lempy.LocalAnalysis(trajectory)
la.print_to_file()
```

With the default options these commands will generate the following files:
* `lem_I_local.dat` and `lem_J_local.dat` contain the instantaneous shape and volume fluctuations for each of the points stored in the trajectory, with `lem_I_local_pmf.dat` and `lem_J_local_pmf.dat` storing the associated normalised histograms.
* `lem_I_local_avg.dat` and `lem_J_local_avg.dat` contain the global instantaneous shape and volume fluctuations of each configuration, computed by averaging the local deformation gradient tensor. `lem_I_local_avg_pmf.dat` and `lem_J_local_avg_pmf.dat` store the associated normalised histograms.
	
The `LocalAnalysis` object takes the following optional arguments that can also change the quantities that are computed:
* `id_centre=None`: sets the given particle as the central bead.
* `id_particles=None`: only uses the particles given to compute the local elastic fluctuations.
* `find_correspondence=True`: if `False` the analysis will assume that points in different configurations have a different identity and will attempt to assign the right identity to each point by using an [iterative closest point](https://en.wikipedia.org/wiki/Iterative_closest_point) construction.
* `relative_to_centre=False`: if set to `True` the analysis will also compute the *global* elastic fluctuations of each configuration by using either the central bead (if `id_centre` is not `None`) or the configuration's centre of mass as reference. If set to `True` the code will print four additional files:
	* `lem_I_global.dat` and `lem_J_global.dat` will contain the global shape and volume fluctuations.
	* `lem_I_global_pmf.dat` and `lem_J_global_pmf.dat` the associated normalised histograms.
	
### Output options

Both `ConvexHullAnalysis`'s and `LocalAnalysis`'s `print_to_file` methods accept an optional parameter `prefix` (which defaults to `ch_` and `lem_`, respectively) that will be prefixed to the name of all output files.

## A slightly more complicated example

This is an example that will analyse a LAMMPS trajectory of a star polymers made of 10 arms of 250 monomers each, where we want to compute the elastic properties of the tips of the arms only:  

```Python
import baggianalysis as ba
import lempy

arms = 10
N_per_arm = 250
id_particles = [a * N_per_arm + 251 for a in range(arms)]
id_centre = 1

parser = ba.LAMMPSDumpParser()
trajectory = ba.FullTrajectory(parser)
trajectory.initialise_from_trajectory_file("f_10_n_250.lammpstrj")

ch = lempy.ConvexHullAnalysis(trajectory)
ch.print_to_file()

la = lempy.LocalAnalysis(trajectory, id_centre=id_centre, id_particles=id_particles, find_correspondence=True, relative_to_centre=True)
la.print_to_file()
```

## Requirements

The package requires `numpy` to perform the required matrix algebra and [`baggianalysis`](https://github.com/lorenzo-rovigatti/baggianalysis) to read the simulation trajectories and perform the convex hull construction. 
