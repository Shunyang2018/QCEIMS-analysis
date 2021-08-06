import numpy as np
import warnings


def centroid_spec(spec, ms2_ppm=None, ms2_da=0.05):
    """
    If both ms2_ppm and ms2_da is defined, ms2_da will be used.
    """
    intensity_order = np.argsort(-spec[:, 1])
    spec_new = []

    # Fast check is the spectrum need centroid.
    mz_array = spec[:, 0]
    need_centroid = 0
    if mz_array.shape[0] > 1:
        mz_delta = mz_array[1:] - mz_array[:-1]
        if ms2_da is not None:
            if np.min(mz_delta) <= ms2_da:
                need_centroid = 1
        else:
            if np.min(mz_delta / mz_array[1:] * 1e6) <= ms2_ppm:
                need_centroid = 1

    if need_centroid:
        for i in intensity_order:
            if ms2_da is None:
                if ms2_ppm is None:
                    raise RuntimeError("MS2 tolerance not defined.")
                else:
                    mz_delta_allowed = ms2_ppm * 1e-6 * spec[i, 0]
            else:
                mz_delta_allowed = ms2_da

            if spec[i, 1] > 0:
                # Find left board for current peak
                i_left = i - 1
                while i_left >= 0:
                    mz_delta_left = spec[i, 0] - spec[i_left, 0]
                    if mz_delta_left <= mz_delta_allowed:
                        i_left -= 1
                    else:
                        break
                i_left += 1

                # Find right board for current peak
                i_right = i + 1
                while i_right < spec.shape[0]:
                    mz_delta_right = spec[i_right, 0] - spec[i, 0]
                    if mz_delta_right <= mz_delta_allowed:
                        i_right += 1
                    else:
                        break

                # Merge those peaks
                intensity_sum = np.sum(spec[i_left:i_right, 1])
                intensity_weighted_sum = np.sum(spec[i_left:i_right, 0] * spec[i_left:i_right, 1])

                spec_new.append([intensity_weighted_sum / intensity_sum, intensity_sum])
                spec[i_left:i_right, 1] = 0

        spec_new = np.array(spec_new)
        # Sort by m/z
        spec_new = spec_new[np.argsort(spec_new[:, 0])]
        return spec_new
    else:
        return spec


def normalize_spec(spec: list, ms2_ppm=None, ms2_da=0.05):
    """
    If both ms2_ppm and ms2_da is defined, ms2_da will be used.
    """
    # Remove intensity==0
    spec = [x for x in spec if x[1] > 0]

    # Convert to numpy array
    spec = np.array(spec, dtype=np.float64)

    if spec.shape[0] > 0:
        # Sort by m/z
        spec = spec[np.argsort(spec[:, 0])]

        # Centroid spectrum
        spec = centroid_spec(spec, ms2_ppm, ms2_da)

        # Normalize the spectrum to sum(intensity)==1
        spec_sum = np.sum(spec[:, 1])
        spec[:, 1] /= spec_sum

    return spec


def match_spec(spec_a, spec_b, ms2_ppm=None, ms2_da=None) -> [[]]:
    """
    Match two spectra, find common peaks. If both ms2_ppm and ms2_da is defined, ms2_da will be used.
    :return: list. Each element in the list is a list contain three elements:
                              m/z, intensity from spec 1; intensity from spec 2.
    """
    a = 0
    b = 0

    spec_merged = []
    peak_b_int = 0.

    while a < spec_a.shape[0] and b < spec_b.shape[0]:
        if ms2_da is not None:
            ms2_ppm = ms2_da / spec_a[a, 0] * 1e6
        mass_delta_ppm = (spec_a[a, 0] - spec_b[b, 0]) / spec_a[a, 0] * 1e6

        if mass_delta_ppm < -ms2_ppm:
            # Peak only existed in spec a.
            spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
            peak_b_int = 0.
            a += 1
        elif mass_delta_ppm > ms2_ppm:
            # Peak only existed in spec b.
            spec_merged.append([spec_b[b, 0], 0., spec_b[b, 1]])
            b += 1
        else:
            # Peak existed in both spec.
            peak_b_int += spec_b[b, 1]
            b += 1

    if peak_b_int > 0.:
        spec_merged.append([spec_a[a, 0], spec_a[a, 1], peak_b_int])
        peak_b_int = 0.
        a += 1

    if b < spec_b.shape[0]:
        spec_merged += [[x[0], 0., x[1]] for x in spec_b[b:]]

    if a < spec_a.shape[0]:
        spec_merged += [[x[0], x[1], 0.] for x in spec_a[a:]]

    if spec_merged:
        spec_merged = np.array(spec_merged, dtype=np.float64)
    else:
        spec_merged = np.array([[0., 0., 0.]], dtype=np.float64)
    return spec_merged


def normalize_distance(dist, dist_range):
    if dist_range[1] == np.inf:
        if dist_range[0] == 0:
            result = 1 - 1 / (1 + dist)
        elif dist_range[1] == 1:
            result = 1 - 1 / dist
        else:
            raise NotImplementedError()
    elif dist_range[0] == -np.inf:
        if dist_range[1] == 0:
            result = -1 / (-1 + dist)
        else:
            raise NotImplementedError()
    else:
        result = (dist - dist_range[0]) / (dist_range[1] - dist_range[0])

    if result < 0:
        result = 0.
    elif result > 1:
        result = 1.

    return result
