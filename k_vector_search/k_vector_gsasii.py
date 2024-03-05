from collections import defaultdict
import GSASIIlattice as G2lat
import random as rand
import sys
import numpy as np

def unique_id_gen(string_list: list) -> list:
    """Generate unique IDs for strings included in the string list and the same
    string will be assigned with the same ID.

    :param string_list: list of strings
    :return: list of integer IDs
    """
    unique_integer = 1
    replaced_strings = {}
    output_list = []

    for string in string_list:
        if string in replaced_strings:
            output_list.append(replaced_strings[string])
        else:
            replaced_strings[string] = unique_integer
            output_list.append(unique_integer)
            unique_integer += 1

    return output_list

def lat_params_to_vec(lat_params: list) -> list:
    """Construct lattice vectors from lattice parameters, according to the
    convention as detailed in the following post,

    https://iris2020.net/2024-03-04-latt_params_to_latt_vecs/

    :param lat_params: list of lattice parameters a, b, c, alpha, beta
                       and gamma
    :return: lattice vectors in the list form, namely, the a, b and c lattice
             vectors given in the Cartesian coordinate
    """
    c = [0., 0., lat_params[2]]
    b = [
        0.,
        lat_params[1] * np.sin(lat_params[3]),
        lat_params[1] * np.cos(lat_params[3])
    ]
    az = lat_params[0] * np.cos(lat_params[4])
    cos_al = np.cos(lat_params[3])
    cos_be = np.cos(lat_params[4])
    cos_ga = np.cos(lat_params[5])
    sin_al = np.sin(lat_params[3])
    top = lat_params[0] * (cos_ga - cos_al * cos_be)
    bottom = sin_al
    ay = top / bottom
    ax = np.sqrt(lat_params[0]**2. - ay**2. - az**2.)
    a = [ax, zy, za]

    return [a, b, c]

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
    nuc_peaks.append(list(ref[:3]) + [ref[4]])

# worry about the phase
# |||||
# vvvvv
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
lat_type = Phase["General"]["SGData"]["SGLatt"].upper()
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

atom_ids = unique_id_gen(atom_types)

# grab the parent unit cell and construct the lattice vectors
cell_params = newPhase["General"]["Cell"][1:7]
lat_vectors = lat_params_to_vec(cell_params)

# TODO: need interface for users to specify the threshold, i.e., delta_d/d
threshold = 0.008

# TODO: Feed the following inputs into the k vector search routine,
#
# lat_vectors
# atom_coords
# atom_ids
# nuc_peaks
# xtra_peaks_d
# threshold
#


