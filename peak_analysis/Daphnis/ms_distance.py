import numpy as np
from .tools import normalize_spec, match_spec, normalize_distance
import scipy.stats


def _match_spec_with_mz_info(spec_a, spec_b, ms2_ppm=None, ms2_da=None):
    """
    Match two spectra, find common peaks. If both ms2_ppm and ms2_da is defined, ms2_da will be used.
    :return: list. Each element in the list is a list contain three elements:
                              m/z from spec 1; intensity from spec 1; m/z from spec 2; intensity from spec 2.
    """
    a = 0
    b = 0

    spec_merged = []
    peak_b_mz = 0.
    peak_b_int = 0.

    while a < spec_a.shape[0] and b < spec_b.shape[0]:
        mass_delta_ppm = (spec_a[a, 0] - spec_b[b, 0]) / spec_a[a, 0] * 1e6
        if ms2_da is not None:
            ms2_ppm = ms2_da / spec_a[a, 0] * 1e6
        if mass_delta_ppm < -ms2_ppm:
            # Peak only existed in spec a.
            spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_mz, peak_b_int])
            peak_b_mz = 0.
            peak_b_int = 0.
            a += 1
        elif mass_delta_ppm > ms2_ppm:
            # Peak only existed in spec b.
            spec_merged.append([0., 0., spec_b[b, 0], spec_b[b, 1]])
            b += 1
        else:
            # Peak existed in both spec.
            peak_b_mz = ((peak_b_mz * peak_b_int) + (spec_b[b, 0] * spec_b[b, 1])) / (peak_b_int + spec_b[b, 1])
            peak_b_int += spec_b[b, 1]
            b += 1

    if peak_b_int > 0.:
        spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_mz, peak_b_int])
        peak_b_mz = 0.
        peak_b_int = 0.
        a += 1

    if b < spec_b.shape[0]:
        spec_merged += [[0., 0., x[0], x[1]] for x in spec_b[b:]]

    if a < spec_a.shape[0]:
        spec_merged += [[x[0], x[1], 0., 0.] for x in spec_a[a:]]

    if spec_merged:
        spec_merged = np.array(spec_merged, dtype=np.float64)
    else:
        spec_merged = np.array([[0., 0., 0., 0.]], dtype=np.float64)
    return spec_merged


def ms_for_id_v1_distance(spec_match):
    r"""
    MSforID distance version 1:

    .. math::

        Similarity = \frac{N_m^4}{N_qN_r(\sum|I_{q,i}-I_{r,i}|)^a}\ ,\ a=0.25

        Distance = \frac{1}{Similarity}

    :math:`N_m`: number of matching fragments, :math:`N_q, N_r`: number of fragments for spectrum p,q
    :return: :math:`Distance`
    """

    i_q = spec_match[:, 1]
    i_r = spec_match[:, 2]

    n_m = np.sum(np.bitwise_and(i_q > 0, i_r > 0))
    n_q = np.sum(i_q > 0)
    n_r = np.sum(i_r > 0)

    a = 0.25
    x = n_m ** 4
    y = n_q * n_r * np.power(np.sum(np.abs(i_q - i_r)), a)

    if x == 0:
        dist = np.inf
    else:
        dist = y / x

    return dist


def ms_for_id_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    r"""
    MSforID distance:

    .. math::

        -\frac{N_m^b(\sum I_{q,i}+2\sum I_{r,i})^c}{(N_q+2N_r)^d+\sum|I_{q,i}-I_{r,i}|+\sum|M_{q,i}-M_{r,i}|},\ \ b=4,\ c=1.25,\ d=2

    The peaks have been filtered with intensity > 0.05.

    :math:`N_m`: number of matching fragments,

    :math:`N_q, N_r`: number of fragments for spectrum p,q,

    :math:`M_q,M_r`: m/z of peak in query and reference spectrum,

    :math:`I_q,I_r`: intensity of peak in query and reference spectrum
    """

    # Filter spectrum to have intensity >0.05
    if len(spec_query) == 0 or len(spec_reference) == 0:
        return np.inf

    spec_query = spec_query[spec_query[:, 1] > 0.05]
    spec_reference = spec_reference[spec_reference[:, 1] > 0.05]

    spec_matched = _match_spec_with_mz_info(spec_query, spec_reference, ms2_ppm, ms2_da)

    b = 4
    c = 1.25
    d = 2

    i_q = spec_matched[:, 1]
    i_r = spec_matched[:, 3]
    matched_peak = np.bitwise_and(i_q > 0, i_r > 0)
    n_m = np.sum(matched_peak)
    n_q = np.sum(i_q > 0)
    n_r = np.sum(i_r > 0)
    i_delta = (i_q - i_r)[matched_peak]
    m_delta = (spec_matched[:, 0] - spec_matched[:, 2])[matched_peak]

    s1 = np.power(n_m, b) * np.power(np.sum(i_q) + 2 * np.sum(i_r), c)
    s2 = np.power(n_q + 2 * n_r, d) + \
         np.sum(np.abs(i_delta)) + \
         np.sum(np.abs(m_delta))

    if s2 == 0:
        similarity = 0.
    else:
        similarity = s1 / s2
    return -similarity


def nist_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    spec_matched = _match_spec_with_mz_info(spec_query, spec_reference, ms2_ppm, ms2_da)

    m_q = spec_matched[:, 0]
    i_q = spec_matched[:, 1]
    m_r = spec_matched[:, 2]
    i_r = spec_matched[:, 3]

    matched_peak = np.bitwise_and(i_q > 0, i_r > 0)
    n_m = np.sum(matched_peak)
    n_q = np.sum(i_q > 0)
    n_r = np.sum(i_r > 0)

    k = 0.6,
    l = 3
    w_q = np.power(i_q, k) * np.power(m_q, l)
    w_r = np.power(i_r, k) * np.power(m_r, l)
    f_d = np.power(np.sum(w_q * w_r), 2) / (np.sum(np.power(w_q, 2)) * np.sum(np.power(w_r, 2)))

    matched_peak_id = np.array(range(matched_peak.shape[0]))[matched_peak]
    f_r = 0
    for i in matched_peak_id[1:]:
        if i_q[i - 1] > 0:
            r_q = i_q[i] / i_q[i - 1]
        else:
            r_q = np.inf

        if i_r[i - 1] > 0:
            r_r = i_r[i] / i_r[i - 1]
        else:
            r_r = np.inf

        if r_q > r_r:
            r = r_r / r_q
        else:
            r = r_q / r_r

        f_r += r

    similarity = (n_q * f_d + f_r) / (n_q + n_m)
    return 1 - similarity


def weighted_dot_product_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    r"""
    Weighted Dot-Product distance:

    .. math::
        TODO:
        1 - \frac{(\sum{Q_iP_i})^2}{\sum{Q_i^2\sum P_i^2}}
    """

    spec_matched = _match_spec_with_mz_info(spec_query, spec_reference, ms2_ppm, ms2_da)
    m_q = spec_matched[:, 0]
    i_q = spec_matched[:, 1]
    m_r = spec_matched[:, 2]
    i_r = spec_matched[:, 3]
    k = 0.6,
    l = 3
    w_q = np.power(i_q, k) * np.power(m_q, l)
    w_r = np.power(i_r, k) * np.power(m_r, l)

    from . import math_distance
    return math_distance.dot_product_distance(w_q, w_r)


def msdial_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    from distance_methods.msdial_helper.pyspec.parser.pymzl.msms_spectrum import MSMSSpectrum

    # Make a fake precursor, as we do not consider the precursor mz.
    precursor_mz = max(np.max(spec_query[:, 0]), np.max(spec_reference[:, 0])) + 1

    if ms2_da is None:
        ms2_da = (precursor_mz - 1) * ms2_ppm * 1e-6

    spec_query = spec_query.tolist()
    spec_reference = spec_reference.tolist()
    a = MSMSSpectrum([(x[0], x[1]) for x in spec_query], precursor_mz=precursor_mz)
    b = MSMSSpectrum([(x[0], x[1]) for x in spec_reference], precursor_mz=precursor_mz)
    return 1 - a.total_similarity(b, 0.1, ms2_da)[-1]


def entropy2_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    spec_merge = np.concatenate([spec_query, spec_reference])
    spec_merge = normalize_spec(spec_merge)
    entropy_increase = 2 * scipy.stats.entropy(spec_merge[:, 1]) - \
                       scipy.stats.entropy(spec_query[:, 1]) - scipy.stats.entropy(spec_reference[:, 1])
    return entropy_increase


def entropy3_distance(spec_query, spec_reference, ms2_ppm=None, ms2_da=None):
    def norm_spec(spec):
        max_intensity = np.max(spec[:, 1])
        # peaks smaller than 0.5% BPI are not useful for most similarity comparisons
        spec = np.array([x for x in spec if x[1] > 0.005 * max_intensity])
        spec_sum = np.sum(spec[:, 1])
        spec[:, 1] /= spec_sum
        return spec

    spec_query = norm_spec(spec_query)
    spec_reference = norm_spec(spec_reference)

    spec_merge = np.concatenate([spec_query, spec_reference])
    spec_merge = normalize_spec(spec_merge)
    entropy_increase = 2 * scipy.stats.entropy(spec_merge[:, 1]) - \
                       scipy.stats.entropy(spec_query[:, 1]) - scipy.stats.entropy(spec_reference[:, 1])
    return entropy_increase
