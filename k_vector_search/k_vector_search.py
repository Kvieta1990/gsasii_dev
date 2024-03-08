#
# ------------------
# k_vector_search.py
# ------------------
#
# This Python file contains the definition of the `kVector` class for the
# identification of optimal k vector(s) for explaining the experimentally
# observed satellite peaks in powder diffraction patterns. Provided the
# refined nucleus structure and the extracted positions of those
# satellite peaks, we first identify the optimal k-path to search along.
# For such a purpose, we were using the k-path finder as reported by
# Y. Hinuma, et al.,
#
# -------------------------------------------------
# http://dx.doi.org/10.1016/j.commatsci.2016.10.015
# -------------------------------------------------
#
# We first search those high symmetry points along the suggested
# k-path, then along the path and finally, if specified to, search
# across the whole first Brillouin zone.
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# Yuanpeng Zhang & Joe Paddison @ Mar-03-2024
# SNS-HFIR, ORNL
# -------------------------------------------
# Vasile Garlea and Stuart Calder are
# acknowledged for their useful comments.
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
import seekpath
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
                       and gamma. Angles should be given in degree.
    :return: lattice vectors in the list form, namely, the a, b and c lattice
             vectors given in the Cartesian coordinate
    """
    a_norm = lat_params[0]
    b_norm = lat_params[1]
    c_norm = lat_params[2]
    alpha = np.deg2rad(lat_params[3])
    beta = np.deg2rad(lat_params[4])
    gamma = np.deg2rad(lat_params[5])

    c = [0., 0., c_norm]
    b = [
        0.,
        b_norm * np.sin(alpha),
        b_norm * np.cos(alpha)
    ]
    az = a_norm * np.cos(beta)
    cos_al = np.cos(alpha)
    cos_be = np.cos(beta)
    cos_ga = np.cos(gamma)
    sin_al = np.sin(alpha)
    top = a_norm * (cos_ga - cos_al * cos_be)
    bottom = sin_al
    ay = top / bottom
    ax = np.sqrt(a_norm**2. - ay**2. - az**2.)
    a = [ax, ay, az]

    return [a, b, c]


class kVector:
    """For k-vector search, given the input structure, the nucleus and
    satellite diffraction peaks.

    :param bravfSym: Bravais lattice symbol
    :param cell: a `3x3` list of floats (cell[0] is the first lattice vector,
                 etc.) The cell vectors should be corresponding to the
                 primitive following the ITA convention.
    :param positions: a `Nx3` list of floats with the atomic coordinates in
                      fractional coordinates (i.e., w.r.t. the cell vectors)
    :param numbers: a length-`N` list with integers identifying uniquely the
                    atoms in the cell
    :param nucPeaks: a `mx4` list for holding the nucleus diffraction peaks
                     (nucPeaks[0] has 4 entries, giving the hkl indeces
                     corresponding to the conventional unit cell setting (ITA),
                     together with the d-spacing)
    :param superPeaks: a length-`n` list for the satellite peak positions in
                       `d`.
    :param threshold: specify the delta_d/d threshold. When the distance
                      between two peaks is smaller than such a threshold, they
                      would then be considered as identical. When searching
                      for the single k vector, if the maximum value among those
                      distances between the nominal and observed positions of
                      the satellite peaks is smaller than the threshold, the
                      corresponding k vector will be returned directly as the
                      optimal solution.
    :param option: (optional) control the scope for k vector search in the
                   Brillouin zone, `0` for high symmetry points only, `1` for
                   high symmetry and edges, `2` for the whole Brillouin zone.
                   Default: 0
    :param kstep: (optional) step of k values for searching along the k path
                  or across the whole 1st Brillouin zone.
                  Default: 0.01
    """
    transMatrix = {
        "P": np.array(
            [
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1]
            ]
        ),
        "cF": 1 / 2 * np.array(
            [
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 0]
            ]
        ),
        "oF": 1 / 2 * np.array(
            [
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 0]
            ]
        ),
        "cI": 1 / 2 * np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1]
            ]
        ),
        "tI": 1 / 2 * np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1]
            ]
        ),
        "oI": 1 / 2 * np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1]
            ]
        ),
        "hR": 1 / 3 * np.array(
            [
                [2, -1, -1],
                [1, 1, -2],
                [1, 1, 1]
            ]
        ),
        "oC": 1 / 2 * np.array(
            [
                [1, 1, 0],
                [-1, 1, 0],
                [0, 0, 2]
            ]
        ),
        "oA": 1 / 2 * np.array(
            [
                [0, 0, 2],
                [1, 1, 0],
                [-1, 1, 0]
            ]
        ),
        "mC": 1 / 2 * np.array(
            [
                [1, -1, 0],
                [1, -1, 1],
                [1, 1, -1]
            ]
        )
    }

    def __init__(self, bravfSym: str, cell: list, positions: list,
                 numbers: list, nucPeaks: list, superPeaks: list,
                 threshold: float, option: int = 0, kstep: float = 0.01):
        self.bravfSym = bravfSym
        self.cell = cell
        self.positions = positions
        self.numbers = numbers
        self.nucPeaks = nucPeaks
        self.superPeaks = superPeaks
        self.threshold = threshold
        self.option = option
        self.kstep = kstep

    def kpathFinder(self) -> dict:
        """Provided the structure inputs, the routine will be collecting the
        inputs into a structure tuple which will be fed into the `get_path`
        routine in the `seekpath` module. The special k points and the k-path
        will be returned.

        :return: the dictionary containing the k-points, k-path and the
                 reciprocal space primitive lattice vectors.
        """
        structure = (self.cell, self.positions, self.numbers)

        # the inputs for the `seekpath.get_path` routine assumes no symmetry
        # and internally in the routine, the symmetry would be identified
        # automatically using the `spglib` module. For the symmetry
        # identification, we need to specify the tolerance (for atomic
        # coordinates, etc.). Here, for our purpose, we know exactly what
        # the symmetry is for our input structure, we were to tune the
        # tolerance value until the expected Bravais lattice type is
        # identified.
        sym_tol = 1.E-5
        while True:
            k_info_tmp = seekpath.get_path(structure, True, symprec=sym_tol)
            if k_info_tmp["bravais_lattice"] == self.bravfSym:
                break
            else:
                sym_tol *= 2.

        k_info = {
            "point_coords": k_info_tmp["point_coords"],
            "path": k_info_tmp["path"],
            "reciprocal_primitive_lattice": k_info_tmp[
                "reciprocal_primitive_lattice"
            ]
        }

        return k_info

    def hklConvToPrim(self, hkl: list) -> list:
        """Convert the hkl indeces in the conventional cell setting to the
        primitive cell setting.

        :param hkl: input hkl indeces in the conventional cell setting
        :return: a list containing the hkl indeces in the primitive cell
                 setting.
        """
        prim_hkl = np.matmul(
            np.array([hkl]),
            kVector.transMatrix[self.bravfSym]
        )

        return list(prim_hkl[0])

    def kVecPrimToConv(self, k_vec: list) -> list:
        """Convert the k vector in the reciprocal primitive lattice setting to
        that in the conventional cell setting.

        :param k_vec: the k vector in the reciprocal primitive lattice setting
        :return: the k vector in the conventional cell setting
        """
        inv_trans_matrix = np.linalg.inv(
            kVector.transMatrix[self.bravfSym]
        )
        k_vec_conv = np.matmul(
            np.array([k_vec]),
            inv_trans_matrix
        )

        return list(k_vec_conv[0])

    def pointOnVector(self, s_point: list, e_point: list,
                      distance: float) -> list:
        """Grab the coordinate of a point on a vector specified by the starting
        and ending points. The distance from the point on the vector to the
        starting point should be given as the parameter.

        :param s_point: a list for the coordinate of the starting point
        :param e_point: a list for the coordinate of the ending point
        :param distance: the distance away from the starting point
        :return: the coordinate of the point on the vector
        """
        vec_length = np.sqrt(sum((x - y)**2 for x, y in zip(s_point, e_point)))

        xp = s_point[0] + distance / vec_length * (e_point[0] - s_point[0])
        yp = s_point[1] + distance / vec_length * (e_point[1] - s_point[1])
        zp = s_point[2] + distance / vec_length * (e_point[2] - s_point[2])

        return [xp, yp, zp]

    def insIntoSortedList(self, lst: list, new_val: float) -> tuple:
        """Insert a new entry into the ascendingly sorted list.

        :param lst: an ascendingly sorted list
        :param new_val: the new entry to be inserted into the sorted list
        :return: tuple containing the new sorted list and the index of the new
                entry in the new sorted list.
        """
        ind = 0
        for i in range(len(lst)):
            if lst[i] < new_val:
                ind = i + 1
            else:
                break
        lst.insert(ind, new_val)

        return (lst, ind)

    def updateCandidateList(self, kpoint: list, k_opt_list: list,
                            k_opt_dist: list, try_neg: bool) -> tuple:
        """For a given k point, we want to first cycle through all the
        satellite peaks. For each satellite peak, we want to find out which
        nucleus peak that the satellite peak should be attached to, meaning the
        nucleus peak for which the hkl plus the trial k vector would yield the
        minimum distance from the observed peak position of the satellite peak
        under consideration. The maximum value among the minimum distances for
        all the satellite peaks will be taken as the indicator for whether the
        given k point is a good candidate. The new `indicator distance` will be
        inserted into the indicator distance list of those top candidates and
        the given k point will be inserted into the corresponding location in
        the top candidates list of the k vectors. If the `indicator distance`
        happens to be smaller than the uncertainty of peak positions (which is
        determined by the instrument resolution), only the top candidate will
        be returned.

        :param kpoint: the trial k vector
        :param k_opt_list: list of top candidates of k vectors
        :param k_opt_dist: list of the `indicator distances` of those top
                           candidates of k vectors
        :param try_neg: flag for controlling whether to try the negative k
                        vector corresponding to the input k vector trial. For
                        the corners of the irreducible Brillouin zone wedge, we
                        don't need to try their negatives, as the k vector, in
                        that case, is different from its negative by a lattice
                        vector in reciprocal space and thus they are equivalent
                        to each other.
        :return: tuple of the updated list of top candidates of k vectors and
                 the updated list of the `indicator distances`
        """
        rep_prim_latt = self.kpathFinder()["reciprocal_primitive_lattice"]

        speaks_closest_dist = list()
        speaks_closest_hkl = list()
        for sp in self.superPeaks:
            dist_min = np.Inf
            hkl_min = list()
            for nucp in self.nucPeaks:
                hkl_conv = nucp[:3]
                hkl_prim = self.hklConvToPrim(hkl_conv)
                sp_norm = np.array(
                    [
                        [
                            hkl_prim[0] + kpoint[0],
                            hkl_prim[1] + kpoint[1],
                            hkl_prim[2] + kpoint[2]
                        ]
                    ]
                )
                sp_norm_cart = np.matmul(
                    sp_norm,
                    rep_prim_latt
                )
                sp_norm_d = 2. * np.pi / np.linalg.norm(sp_norm_cart)
                sp_rel_diff = abs(sp_norm_d - sp) / sp
                if sp_rel_diff < dist_min:
                    dist_min = sp_rel_diff
                    hkl_min = hkl_conv
                if try_neg:
                    kpoint_neg = [-item for item in kpoint]
                    sp_norm = np.array(
                        [
                            [
                                hkl_prim[0] + kpoint_neg[0],
                                hkl_prim[1] + kpoint_neg[1],
                                hkl_prim[2] + kpoint_neg[2]
                            ]
                        ]
                    )
                    sp_norm_cart = np.matmul(
                        sp_norm,
                        rep_prim_latt
                    )
                    sp_norm_d = 2. * np.pi / np.linalg.norm(sp_norm_cart)
                    sp_rel_diff = abs(sp_norm_d - sp) / sp
                    if sp_rel_diff < dist_min:
                        dist_min = sp_rel_diff
                        hkl_min = hkl_conv
            speaks_closest_dist.append(dist_min)
            speaks_closest_hkl.append(hkl_min)

        if max(speaks_closest_dist) <= self.threshold:
            k_opt_list = [kpoint]
            k_opt_dist = [max(speaks_closest_dist)]
            return (k_opt_list, k_opt_dist)
        else:
            k_opt_new = self.insIntoSortedList(
                k_opt_dist,
                max(speaks_closest_dist)
            )
            k_opt_list.insert(k_opt_new[1], kpoint)

            return (k_opt_list[:10], k_opt_new[0][:10])

    def kOptFinder(self) -> list:
        """This is the kernel of the class, defining the method for searching
        over the Brillouin zone for optimal k vector that best explains the
        observed satellite peaks, given the already refined nucleus structure.

        :return: the list of top candidates of the optimal k vector (in the
                 reciprocal `primitive cell` setting). If the top one is below
                 the threshold, i.e., the distance between the nominal and
                 observed positions of those satellite peaks smaller than the
                 uncertainty of peak positions (which is determined by the
                 instrument resolution), only one candiate will be returned.
                 Otherwise, top 10 candidates will be returned.
        """
        hs_points = self.kpathFinder()["point_coords"]
        rep_prim_latt = self.kpathFinder()["reciprocal_primitive_lattice"]

        k_opt_list = list()
        k_opt_dist = list()

        if len(self.nucPeaks) < np.floor(len(self.superPeaks) / 2):
            err_msg = "The number of nucleus peaks is not enough for "
            err_msg += "k vector determination."
            raise ValueError("err_msg")
        else:
            if len(self.nucPeaks) >= len(self.superPeaks):
                # search over those high symmetry points on the suggested k
                # path In this case, we don't need to consider their negatives
                # as they are equivalent.
                for _, kpoint in hs_points.items():
                    k_opt_tmp = self.updateCandidateList(
                        kpoint,
                        k_opt_list,
                        k_opt_dist,
                        False
                    )
                    k_opt_list = k_opt_tmp[0]
                    k_opt_dist = k_opt_tmp[1]

                    if len(k_opt_list) == 1:
                        return k_opt_list

                if self.option == 1:
                    # search along the k-path
                    k_paths = self.kpathFinder()
                    for k_path in k_paths["path"]:
                        seg_len = self.kstep
                        k_path_s = k_paths["point_coords"][k_path[0]]
                        k_path_e = k_paths["point_coords"][k_path[1]]
                        k_path_vec = np.array(
                            [
                                [y - x] for y, x in zip(k_path_e, k_path_s)
                            ]
                        )
                        k_path_vec_cart = np.matmul(
                            k_path_vec,
                            rep_prim_latt
                        )
                        k_path_len = np.linalg.norm(k_path_vec_cart)
                        while seg_len < k_path_len:
                            kpoint = self.pointOnVector(
                                k_path_s,
                                k_path_e,
                                seg_len
                            )
                            k_opt_tmp = self.updateCandidateList(
                                kpoint,
                                k_opt_list,
                                k_opt_dist,
                                True
                            )
                            k_opt_list = k_opt_tmp[0]
                            k_opt_dist = k_opt_tmp[1]

                            if len(k_opt_list) == 1:
                                return k_opt_list
                elif self.option == 2:
                    # search over the whole 1st Brillouin zone
                    grid_num = int(0.5 / self.kstep) - 1
                    for i in range(grid_num):
                        for j in range(grid_num):
                            for m in range(grid_num):
                                h_tmp = self.kstep * (i + 1)
                                k_tmp = self.kstep * (j + 1)
                                l_tmp = self.kstep * (m + 1)
                                kpoint = [h_tmp, k_tmp, l_tmp]
                                k_opt_tmp = self.updateCandidateList(
                                    kpoint,
                                    k_opt_list,
                                    k_opt_dist,
                                    True
                                )
                                k_opt_list = k_opt_tmp[0]
                                k_opt_dist = k_opt_tmp[1]

                                if len(k_opt_list) == 1:
                                    return k_opt_list
                else:
                    # do nothing if specified to search only those high
                    # symmetry points
                    pass

                return k_opt_list
            else:
                if self.option == 0:
                    err_msg = "The number of nucleus peaks is less than that "
                    err_msg += "of the satellite peaks, and thus the search "
                    err_msg += "should go beyond the high symmetry point. \n"
                    err_msg += "Please select to search either along k paths "
                    err_msg += "or over the whole 1st Brillouin zone."
                    raise ValueError(err_msg)
                else:
                    warn_msg = "The number of nucleus peaks is less than that "
                    warn_msg += "of the satellite peaks, thus skipping the "
                    warn_msg += "search over the high symmetry points."
                    print(f"[Warning] {warn_msg}")
