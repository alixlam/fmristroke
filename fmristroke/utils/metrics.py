import itertools

import nibabel as nib
import numpy as np
import pandas as pd
from bids import layout
from nilearn.maskers import NiftiLabelsMasker
from numpy import linalg as LA
from scipy.stats import pearsonr, ranksums


def compute_isc(time_series_sub, axis=0):
    """Compute Inter subject correlation"""
    results = []
    for i, ts1 in enumerate(time_series_sub):
        for j, ts2 in enumerate(time_series_sub):
            if j == i:
                continue
            results.append(_pearsonr(ts1, ts2, axis=axis))
    results = np.array(results)

    # Mean of correlation not possible - arctan first
    result_per_node = np.tanh(np.mean(np.arctanh(results), axis=0))
    return result_per_node


def compute_intersubjectconnectivity(conn_sub, axis=0):
    """Compute Inter subject connectivity"""
    results = []
    for pair in itertools.combinations(conn_sub, r=2):
        results.append(geodesic(*pair))
    results = np.array(results)
    return results


def _ensure_symmetric(Q):
    """
    computation is sometimes not precise (round errors),
    so ensure matrices that are supposed to be
    symmetric are symmetric
    """
    return (Q + np.transpose(Q)) / 2


def _vectorize(Q):
    """
    given a symmetric matrix (FC), return unique
    elements as an array. Ignore diagonal elements
    """
    # extract lower triangular matrix
    tri = np.tril(Q, -1)

    vec = []
    for ii in range(1, tri.shape[0]):
        for jj in range(ii):
            vec.append(tri[ii, jj])

    return np.asarray(vec)


def geodesic(FC1, FC2):
    """
    dist = sqrt(trace(log^2(M)))
    M = Q_1^{-1/2}*Q_2*Q_1^{-1/2}
    """
    FC1 = _ensure_symmetric(FC1)
    FC2 = _ensure_symmetric(FC2)

    # compute Q_1^{-1/2} via eigen value decmposition
    u, s, _ = LA.svd(FC1, full_matrices=True)

    # lift very small eigen values
    for ii, s_ii in enumerate(s):
        if s_ii < 10 ** (-3):
            s[ii] = 10 ** (-3)

    """
        since FC1 is in S+, u = v, u^{-1} = u'
        FC1 = usu^(-1)
        FC1^{1/2} = u[s^{1/2}]u'
        FC1^{-1/2} = u[s^{-1/2}]u'
        """
    FC1_mod = u @ np.diag(s ** (-1 / 2)) @ np.transpose(u)
    M = FC1_mod @ FC2 @ FC1_mod

    """
        trace = sum of eigenvalues;
        np.logm might have round errors,
        implement using svd instead
        """
    _, s, _ = LA.svd(M, full_matrices=True)

    return np.sqrt(np.sum(np.log(s) ** 2))


def compute_run_similarity(ts1, ts2):
    """
    Compute pearson correlation between two runs from same subject
    """
    results_per_node = _pearsonr(ts1, ts2, axis=axis)

    return results_per_node


def compute_lesion_correlation(lesion_time_series, brain_time_series):
    """
    Compute lesion correlation
    Correlation of lesion signal with other voxels
    """

    lesion_to_voxel_correlations = (
        np.dot(brain_time_series.T, lesion_time_series)
        / lesion_time_series.shape[0]
    )

    mean = np.tanh(np.mean(np.arctanh(lesion_to_voxel_correlations)))
    return lesion_to_voxel_correlations, mean


def _pearsonr(x, y, axis=0):
    """
    Pearson correlation between x and y  on chosen axis
    Examples
    --------
    >>> X = np.array([[1, 2, 3, 4, 5, 6, 7] for _ in range(2)])
    >>> Y = np.array([[10, 9, 2.5, 6, 4, 3, 2] for _ in range(2)])
    >>> r = _pearsonr(X,Y, axis=1)
    >>> r.shape
    (2,)
    >>> r[0]
    -0.828503883588428
    """
    # Compute statistics
    x_mean = np.mean(x, axis=axis, keepdims=True)
    y_mean = np.mean(y, axis=axis, keepdims=True)

    x_std = np.sum(np.square(x - x_mean), axis=axis, keepdims=True)
    y_std = np.sum(np.square(y - y_mean), axis=axis, keepdims=True)

    denominator = np.sqrt(x_std * y_std)
    numerator = np.sum(
        (y - y_mean) * (x - x_mean),
        axis=axis,
        keepdims=True,
    )
    return np.mean(np.divide(numerator, denominator), axis=axis)


def compute_FCC(connectivity_mat, atlas_labels):
    """FCC measure from `A multi-measure approach for assessing the performance of fMRI preprocessing strategies in resting-state functional connectivity`, Kassinopoulos et al."""
    WNE, BNE = _get_WNE_BNE(connectivity_mat, atlas_labels)
    # Wilcoxin rank-sum test
    # H0: WNEs and BNEs samples from continuous distribution with equal medians
    # H1: The distribution underlying WNEs is stochastically greater than the distribution underlying BNEs.
    # FCC defined as the Z-statistic of the Wilcoxon rank.
    z_stat, p_value = ranksums(WNE, BNE, "greater")
    return z_stat, p_value


def _get_WNE_BNE(connectivity_mat, atlas_labels):
    """
    Get WNE and BNE from connectivity matrices when roi are grouped by networks
    Examples
    ---------
    >>> con_mat = np.array([[1,2,3,4],[2,1,2,3],[3,2,1,2],[4,3,2,1]])
    >>> atlas_labels = pd.DataFrame({'Network': ['a', 'a', 'b', 'c']})
    >>> WNE, BNE =  _get_WNE_BNE(con_mat, atlas_labels)
    >>> WNE
    [2]
    >>> BNE
    [3, 2, 4, 3, 2]
    """
    # Get networks in atlas labels
    if isinstance(atlas_labels, str):
        atlas_labels = pd.read_csv(atlas_labels)

    assert (
        "Network" in atlas_labels.columns
    ), "atlas labels should contain a 'Network' column"

    networks = set(atlas_labels["Network"])
    WNE = []
    BNE = []
    for row, lab_row in zip(
        range(connectivity_mat.shape[0]), atlas_labels["Network"]
    ):
        for col, lab_col in zip(
            range(connectivity_mat.shape[1]), atlas_labels["Network"]
        ):
            if row > col:
                if lab_col == lab_row:
                    WNE.append(connectivity_mat[row, col])
                else:
                    BNE.append(connectivity_mat[row, col])

    return WNE, BNE


def get_group_intersubconn(data_layout):
    """
    Compute quality check values for group
    """

    queries = {
        "datatype": "func",
        "extension": ".npy",
        "desc": "connectivity",
        "suffix": "mat",
    }

    sessions = data_layout.get_sessions()
    tasks = data_layout.get_tasks()
    pipelines = data_layout.get_pipelines()

    results = []
    for ses in tqdm.tqdm(sessions):
        for task in tqdm.tqdm(tasks):
            for pipeline in tqdm.tqdm(pipelines):
                if pipeline not in [
                    "SimpleLesionGS",
                    "SimpleGS",
                    "IClesionGS",
                ]:
                    pass
                queries["session"] = ses
                queries["task"] = task
                queries["pipeline"] = pipeline
                data = data_layout.get(return_type="filename", **queries)
                conn_mat = [np.load(f) for f in data]
                conn_mat = [mat for mat in conn_mat if mat.shape[0] == 400]
                final_conn = np.array(conn_mat)
                isc = compute_intersubjectconnectivity(final_conn)
                results.append(
                    pd.DataFrame(
                        {
                            "pipeline": [pipeline] * len(isc),
                            "session": [ses] * len(isc),
                            "task": [task] * len(isc),
                            "isc": isc,
                        }
                    )
                )

    results = pd.concat(results)

    return results


def get_group_roicon(data_layout):
    queries = {
        "datatype": "func",
        "extension": ".tsv",
        "desc": "connectivity",
        "suffix": "roi",
    }

    sessions = data_layout.get_sessions()
    tasks = data_layout.get_tasks()

    results = pd.DataFrame()

    for ses in sessions:
        queries["session"] = ses
        for task in tasks:
            queries["task"] = task
            datas = data_layout.get(return_type="filename", **queries)
            for data in datas:
                df = pd.read_table(data, sep="\t")
                df["session"] = [ses] * df.shape[0]
                df["task"] = [task] * df.shape[0]
                results = pd.concat((results, df), axis=0)
    return results


def get_group_significantroicon(data_layout, n):
    queries = {
        "datatype": "func",
        "extension": ".tsv",
        "desc": "connectivity",
        "suffix": "roi",
    }

    sessions = data_layout.get_sessions()
    tasks = data_layout.get_tasks()

    results = pd.DataFrame()
    dist = scipy.stats.beta(n / 2 - 1, n / 2 - 1, loc=-1, scale=2)

    for ses in sessions:
        queries["session"] = ses
        for task in tasks:
            queries["task"] = task
            datas = data_layout.get(return_type="filename", **queries)
            for data in datas:
                df = pd.read_table(data, sep="\t")
                new_df = pd.DataFrame()
                for col in df:
                    p = 2 * dist.cdf(-np.abs(df[col]))

                    n_significant = np.sum(p < 0.05)
                    new_df[col] = [n_significant]
                new_df["session"] = [ses]
                new_df["task"] = [task]
                results = pd.concat((results, new_df), axis=0)
    return results


def get_group_FCC(data_layout):
    queries = {
        "datatype": "func",
        "extension": ".tsv",
        "desc": "connectivity",
        "suffix": "FCC",
    }
    sessions = data_layout.get_sessions()
    tasks = data_layout.get_tasks()

    results = pd.DataFrame()

    for ses in sessions:
        queries["session"] = ses
        for task in tasks:
            queries["task"] = task
            datas = data_layout.get(return_type="filename", **queries)
            for data in datas:
                df = pd.read_table(data, sep="\t")
                df["session"] = [ses] * df.shape[0]
                df["task"] = [task] * df.shape[0]
                results = pd.concat((results, df), axis=0)
    return results


def get_group_intrasubconn(data_layout, atlas):
    queries = {
        "datatype": "func",
        "extension": ".nii.gz",
        "desc": "denoised",
        "suffix": "bold",
        "space": "MNI152NLin2009cAsym",
    }

    sessions = data_layout.get_sessions()
    tasks = data_layout.get_tasks()
    pipelines = data_layout.get_pipelines()

    results = []

    masker = NiftiLabelsMasker(atlas)
    connectivity = ConnectivityMeasure(
        kind="correlation",
    )

    for ses in tqdm.tqdm(sessions):
        for task in tqdm.tqdm(tasks):
            for pipeline in tqdm.tqdm(pipelines):
                queries["session"] = ses
                queries["task"] = task
                queries["pipeline"] = pipeline
                data = data_layout.get(return_type="filename", **queries)
                for d in data:
                    func = nib.load(d)
                    ts = masker.fit_transform(d)
                    i = ts.shape[0]
                    ts1 = ts[: i // 2, :]
                    ts2 = ts[i // 2 :, :]
                    conn1 = connectivity.fit_transform([ts1])[0]
                    conn2 = connectivity.fit_transform([ts2])[0]
                    isc = compute_intersubjectconnectivity(
                        np.array([conn1, conn2])
                    )
                    results.append(
                        pd.DataFrame(
                            {
                                "pipeline": [pipeline] * len(isc),
                                "session": [ses] * len(isc),
                                "task": [task] * len(isc),
                                "isc": isc,
                            }
                        )
                    )

    results = pd.concat(results)

    return results