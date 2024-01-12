import pandas as pd
import numpy as np
from scipy.stats import ranksums

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

def compute_icc(connectivity_mat):
    """
    Compute Intra subject correlation
    
    """
    return

def compute_lesion_correlation(lesion_time_series, brain_time_series):
    """
    Compute lesion correlation
    Correlation of lesion signal with other voxels
    """
    
    lesion_to_voxel_correlations = (
        np.dot(brain_time_series.T, lesion_time_series) / lesion_time_series.shape[0]
    )
    
    mean = np.tanh(
        np.mean(np.arctanh(leion_to_voxel_correlations))
    )
    return lesion_to_voxel_correlations, mean

def _pearsonr(x, y, axis=0):
    """
    Pearson correlation between x and y  on chosen axis

    Examples
    --------   
    >>> X = np.array([[1, 2, 3, 4, 5, 6, 7] for _ in range(2)]
    >>> Y = np.array([[10, 9, 2.5, 6, 4, 3, 2] for _ in range(2)])
    >>> r = pearsonr(X,Y, axis=1)
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
    z_stat, p_value = ranksums(WNE, BNE, 'greater')
    return z_stat

def _get_WNE_BNE(connectivity_mat, atlas_labels):
    """
    Get WNE and BNE from connectivity matrices when roi are grouped by networks  
    
    Examples
    ---------
    >>> con_mat = np.array([
        [1,2,3,4],
        [2,1,2,3],
        [3,2,1,2],
        [4,3,2,1]
    ])
    >>> atlas_labels = pd.DataFrame({'Network'}: ['a', 'a', 'b', 'c'])
    >>> WNE, BNE =  _get_WNE_BNE(conn_mat, atlas_labels)
    >>> WNE
    [2]
    >>> BNE
    [3,2,4,3,2]
    """
    # Get networks in atlas labels
    if isinstance(atlas_labels, str):
        atlas_labels = pd.read_csv(atlas_labels)
    
    assert "Network" in atlas_labels.columns, "atlas labels should contain a 'Network' column"
    
    networks = set(atlas_labels["Network"])
    WNE = []
    BNE = []
    for net in networks:
        for row, lab_row in zip(connectivity_mat.shape[0], atlas_labels["Network"]):
            for col, lab_col in zip(connectivity_mat.shape[1], atlas_labels["Network"]):
                if row > col:
                    if lab_col == lab_row:
                        WNE.append(connectivity_mat[row, col])
                    else:
                        BNE.append(connectivity_mat[row, col])
    
    return WNE, BNE