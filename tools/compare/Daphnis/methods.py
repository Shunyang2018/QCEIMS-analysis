import numpy as np
from . import math_distance, ms_distance
from .tools import normalize_spec, match_spec, normalize_distance

import warnings

methods_name = {
    "entropy": "Entropy distance",
    "weighted_entropy": "Dynamic weighted entropy distance",
    "entropy_reverse": "Reverse entropy distance",
    "euclidean": "Euclidean distance",
    "manhattan": "Manhattan distance",
    "chebyshev": "Chebyshev distance",
    "squared_euclidean": "Squared Euclidean distance",
    "fidelity": "Fidelity distance",
    "matusita": "Matusita distance",
    "squared_chord": "Squared-chord distance",
    "bhattacharya_1": "Bhattacharya 1 distance",
    "bhattacharya_2": "Bhattacharya 2 distance",
    "harmonic_mean": "Harmonic mean distance",
    # "pearson_chi_squared": "Pearson χ2 distance",
    # "neyman_chi_squared": "Neyman χ2 distance",
    "probabilistic_symmetric_chi_squared": "Probabilistic symmetric χ2 distance",
    "topsoe": "Topsøe distance",
    # "chernoff": "Chernoff distance",
    "ruzicka": "Ruzicka distance",
    "roberts": "Roberts distance",
    "intersection": "Intersection distance",
    "motyka": "Motyka distance",
    "canberra": "Canberra distance",
    # "kulczynski_1": "Kulczynski 1 distance",
    "baroni_urbani_buser": "Baroni-Urbani-Buser distance",
    "penrose_size": "Penrose size distance",
    "mean_character": "Mean character distance",
    "lorentzian": "Lorentzian distance",
    "penrose_shape": "Penrose shape distance",
    "clark": "Clark distance",
    "hellinger": "Hellinger distance",
    "whittaker_index_of_association": "Whittaker index of association distance",
    "symmetric_chi_squared": "Symmetric χ2 distance",
    "pearson_correlation": "Pearson/Spearman Correlation Coefficient",
    # "similarity_index": "Similarity Index",
    "improved_similarity": "Improved Similarity",
    "absolute_value": "Absolute Value Distance",
    "dot_product": "Dot-Product (cosine)",
    "dot_product_reverse": "Reverse dot-Product (cosine)",
    "spectral_contrast_angle": "Spectral Contrast Angle",
    "wave_hedges": "Wave Hedges distance",
    "cosine": "Cosine distance",
    "jaccard": "Jaccard distance",
    "dice": "Dice distance",
    "inner_product": "Inner Product distance",
    "divergence": "Divergence distance",
    # "additive_symmetric_chi_squared": "Additive symmetric χ2 distance",
    # "jensen_difference": "Jensen difference",
    # "kumar_johnson": "Kumar-Johnson distance",
    "avg_l": "Avg (L1, L∞) distance",
    # "vicis_wave_hadges": "Vicis-Wave Hadges distance",
    # "vicis_symmetric_chi_squared_1": "Vicis-Symmetric χ2 1 distance",
    # "vicis_symmetric_chi_squared_2": "Vicis-Symmetric χ2 2 distance",
    "vicis_symmetric_chi_squared_3": "Vicis-Symmetric χ2 3 distance",
    # "max_symmetric_chi_squared": "Max-Symmetric χ2 distance",
    # "min_symmetric_chi_squared": "Min-Symmetric χ2 distance",
    # "mahalanobis": "Mahalanobis distance",
    "ms_for_id_v1": "MSforID distance version 1",
    "ms_for_id": "MSforID distance",
    "nist": "NIST distance",
    "weighted_dot_product": "NIST weighted dot product distance",
    # "msdial": "MS-DIAL distance",
    "entropy2": "Entropy v2",
    "entropy3": "Entropy v3"
}

methods_range = {
    "entropy": [0, np.log(4)],
    "weighted_entropy": [0, np.log(4)],
    "entropy2": [0, np.log(4)],
    "entropy3": [0, np.log(4)],
    "entropy_reverse": [0, np.log(4)],
    "absolute_value": [0, 2],
    "avg_l": [0, 1.5],
    "bhattacharya_1": [0, np.arccos(0) ** 2],
    "bhattacharya_2": [0, np.inf],
    "canberra": [0, np.inf],
    "clark": [0, np.inf],
    "divergence": [0, np.inf],
    "euclidean": [0, np.sqrt(2)],
    "hellinger": [0, np.inf],
    "improved_similarity": [0, np.inf],
    # "jensen_difference": "*",
    # "kulczynski_1": [0,np.inf],
    "lorentzian": [0, np.inf],  # TODO:?
    "manhattan": [0, 2],
    "matusita": [0, np.sqrt(2)],
    "mean_character": [0, 2],
    "motyka": [-0.5, 0],
    "ms_for_id": [-np.inf, 0],
    "ms_for_id_v1": [0, np.inf],
    # "neyman_chi_squared": [0,np.inf],
    "pearson_correlation": [-1, 1],
    "penrose_shape": [0, np.sqrt(2)],
    "penrose_size": [0, np.inf],
    "probabilistic_symmetric_chi_squared": [0, 1],
    "similarity_index": [0, np.inf],
    "squared_chord": [0, 2],
    "squared_euclidean": [0, 2],
    "symmetric_chi_squared": [0, 0.5 * np.sqrt(2)],
    "topsoe": [0, np.sqrt(2)],
    # "vicis_symmetric_chi_squared_1": [0,np.inf],
    "vicis_symmetric_chi_squared_3": [0, 2],
    "wave_hedges": [0, np.inf],
    "whittaker_index_of_association": [0, np.inf]
}


def distance(spec_query: np.ndarray, spec_reference: np.ndarray, method: str,
             normalize_result: bool = True, spectrum_refined: bool = False,
             ms2_ppm: float = 50, ms2_da: float = None) -> float:
    """
    Calculate the distance between two spectra, find common peaks.
    If both ms2_ppm and ms2_da is defined, ms2_da will be used.
    :param normalize_result: normalize the result into [0,1] or not.
    :param spectrum_refined: Is the spectra be centroid and normalized to sum of intensity =1 or not.
    """
    if not spectrum_refined:
        spec_query = normalize_spec(spec_query, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
        spec_reference = normalize_spec(spec_reference, ms2_ppm=ms2_ppm, ms2_da=ms2_da)

    # warnings.filterwarnings('error')
    if spec_query.shape[0] > 0 and spec_reference.shape[0] > 0:
        # Calculate distance
        function_name = method + "_distance"
        if hasattr(math_distance, function_name):
            f = getattr(math_distance, function_name)
            spec_matched = match_spec(spec_a=spec_query, spec_b=spec_reference, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
            dist = f(spec_matched[:, 1], spec_matched[:, 2])

        elif hasattr(ms_distance, function_name):
            f = getattr(ms_distance, function_name)
            if method in {"ms_for_id", "nist", "weighted_dot_product", "msdial", "entropy2", "entropy3"}:
                dist = f(spec_query, spec_reference, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
            else:
                spec_matched = match_spec(spec_a=spec_query, spec_b=spec_reference, ms2_ppm=ms2_ppm, ms2_da=ms2_da)
                dist = f(spec_matched)
        else:
            raise RuntimeError("Method name error!")

        # Normalize result to [0,1]
        if normalize_result:
            if method not in methods_range:
                dist_range = [0, 1]
            else:
                dist_range = methods_range[method]

            dist = normalize_distance(dist, dist_range)
        return dist

    else:
        if normalize_result:
            return 1
        else:
            return np.inf


def all_distance(spec_query_raw, spec_reference_raw, ms2_ppm=50):
    result = {}
    for method in methods_name:
        dist = distance(spec_query_raw, spec_reference_raw, method=method, ms2_ppm=ms2_ppm)
        result[method] = float(dist)
        # result[methods_name[method]] = float(dist)
    return result
