import k_vector_search as kvs
import GSASIIlattice as G2lat
import random as rand
import sys
import numpy as np
import copy

# grab those extra peaks and obtain the corresponding d-spacing
Id = G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId, 'Instrument Parameters')
Parms, _ = G2frame.GPXtree.GetItemPyData(Id)

xtra_peaks_d = list()
for extra_peak in data["xtraPeaks"]:
    dsp_tmp = G2lat.Pos2dsp(Parms, extra_peak[0])
    xtra_peaks_d.append(dsp_tmp)

# grab the nucleus reflections
#
# TODO: Need user selection for the nucleus phase which the magnetic
# structure is attached to. Here, as the draft code, we take the first
# phase in the list to proceed. In the final codes, we only need to worry
# about the `phase_use` variable below.
refDict = G2frame.GPXtree.GetItemPyData(
    G2gd.GetGPXtreeItemId(G2frame, G2frame.PatternId, 'Reflection Lists')
)
phase_use = list(refDict.keys())[0]

nuc_peaks = list()
for ref in refDict[phase_use]["RefList"]:
    # only grab the hkl index and d-spacing
    nuc_peaks.append(list(ref[:3]) + [ref[4]])

# worry about the phase
#
# first, we grab the phase to use
_, Phases = G2frame.GetUsedHistogramsAndPhasesfromTree()
Phase = Phases[phase_use]

# grab the Bravais lattice type
#
# given the lattice type and lattice system, the Bravais lattice type can be
# determined. see the comparison table in the Wikipedia link below for the
# correspondence,
# https://en.wikipedia.org/wiki/Crystal_system
# here, we need to assign the Bravais lattice with specific names as inputs for
# the `seekpath` routine.
lat_type = Phase["General"]["SGData"]["SGLatt"]
lat_sym = Phase["General"]["SGData"]["SGSys"]
if lat_type == "P":
    brav_sym = "P"
else:
    if lat_sym == "trigonal":
        brav_sym = "hR"
    else:
        brav_sym = lat_sym[0] + lat_type

# grab all atomic coordinates in the P1 symmetry
#
# define some matrix as necessary inputs for generating the P1 structure.
Trans = np.eye(3)
Uvec = np.zeros(3)
Vvec = np.zeros(3)

# expand the structure to P1 symmetry
newPhase = copy.deepcopy(Phase)
newPhase['ranId'] = rand.randint(0, sys.maxsize)
newPhase['General']['SGData'] = G2spc.SpcGroup('P 1')[1]
newPhase, _ = G2lat.TransformPhase(Phase, newPhase, Trans, Uvec, Vvec, False)
atoms_pointer = newPhase['General']['AtomPtrs']

atom_coords = list()
atom_types = list()
for atom in newPhase["Atoms"]:
    coord_tmp = atom[atoms_pointer[0]:atoms_pointer[0] + 3]
    atom_coords.append(coord_tmp)
    type_tmp = atom[atoms_pointer[1]]
    atom_types.append(type_tmp)

# this will turn each of the atom types into a unique integer number, which is
# required by the `seekpath` routine
atom_ids = kvs.unique_id_gen(atom_types)

# grab the parent unit cell and construct the lattice vectors
cell_params = newPhase["General"]["Cell"][1:7]
lat_vectors = kvs.lat_params_to_vec(cell_params)

# TODO: need interface for users to specify the threshold for k vector search.
# refer to the doc of the `kVector`class,
#
# https://yr.iris-home.net/kvectordoc
#
threshold = 0.008

# perform the k vector search given all the collected inputs
k_search = kvs.kVector(
    brav_sym,
    lat_vectors,
    atom_coords,
    atom_ids,
    nuc_peaks,
    xtra_peaks_d,
    threshold
)
k_opt = k_search.kOptFinder()
