import numpy as np
import pandas as pd

from sklearn.utils.validation import (
    check_array,
    check_symmetric,
    check_consistent_length,
)

def dist2(X, C):
    """
    Description: Computes the Euclidean distances between all pairs of data point given

    Usage: dist2(X, C)
    X: A data matrix where each row is a different data point
    C: A data matrix where each row is a different data point. If this matrix is the same as X,
    pairwise distances for all data points in X are computed.

    Return: Returns an N x M matrix where N is the number of rows in X and M is the number of rows in C.

    Author: Dr. Anna Goldenberg, Bo Wang, Aziz Mezlini, Feyyaz Demir
    Python Version Rewrite: Rex Ma

    Examples:
        # Data1 is of size n x d_1, where n is the number of patients, d_1 is the number of genes,
        # Data2 is of size n x d_2, where n is the number of patients, d_2 is the number of methylation
        Dist1 = dist2(Data1, Data1)
        Dist2 = dist2(Data2, Data2)
    """

    ndata = X.shape[0]
    ncentres = C.shape[0]

    sumsqX = np.sum(X * X, axis=1)
    sumsqC = np.sum(C * C, axis=1)

    XC = 2 * (np.matmul(X, np.transpose(C)))

    res = (
        np.transpose(np.reshape(np.tile(sumsqX, ncentres), (ncentres, ndata)))
        + np.reshape(np.tile(sumsqC, ndata), (ndata, ncentres))
        - XC
    )

    return res


def _find_dominate_set_relative(W, K=20):
    """
    Retains `K` strongest edges for each sample in `W`
    Parameters
    ----------
    W : (N, N) array_like
        Input data
    K : (0, N) int, optional
        Number of neighbors to retain. Default: 20
    Returns
    -------
    Wk : (N, N) np.ndarray
        Thresholded version of `W`
    """

    # let's not modify W in place
    Wk = W.copy()

    # determine percentile cutoff that will keep only `K` edges for each sample
    # remove everything below this cutoff
    cutoff = 100 - (100 * (K / len(W)))
    Wk[Wk < np.percentile(Wk, cutoff, axis=1, keepdims=True)] = 0

    # normalize by strength of remaining edges
    Wk = Wk / np.nansum(Wk, axis=1, keepdims=True)

    Ws = Wk + np.transpose(Wk)

    return Ws





def _scaling_normalized_pd(W, ratio):
    """
    Adds `alpha` to the diagonal of pandas dataframe `W`

    Parameters
    ----------
    W : (N, N) array_like
        Similarity array from SNF

    Returns
    -------
    W : (N, N) np.ndarray
        Stable-normalized similiarity array
    """

    # W = W.copy()

    # add `alpha` to the diagonal and symmetrize `W`
    rowSum = np.sum(W, 1) - np.diag(W)
    rowSum[rowSum == 0] = 1

    W = (W / rowSum) * 0.5 * ratio 

    W_np = W.to_numpy(copy=True)
    np.fill_diagonal(W_np, 1-0.5*ratio)
    W = pd.DataFrame(W_np, index=W.index, columns=W.columns) 

    W = check_symmetric(W, raise_warning=False)

    return W


def _stable_normalized_pd(W):
    """
    Adds `alpha` to the diagonal of pandas dataframe `W`

    Parameters
    ----------
    W : (N, N) array_like
        Similarity array from SNF

    Returns
    -------
    W : (N, N) np.ndarray
        Stable-normalized similiarity array
    """

    # add `alpha` to the diagonal and symmetrize `W`
    # print(W.isna().any())
    
    rowSum = np.sum(W, 1) - np.diag(W)
    
    # sys.exit()
    rowSum[rowSum == 0] = 1

    W = W / (2 * rowSum)
    # print(W)
   
    W_np = W.to_numpy(copy=True)
    # print(W_np)
    
    np.fill_diagonal(W_np, 0.5)
    W = pd.DataFrame(W_np, index=W.index, columns=W.columns)
    

    W = check_symmetric(W, raise_warning=False)
    
    return W
def _stable_normalized(W):
    """
    Adds `alpha` to the diagonal of `W`

    Parameters
    ----------
    W : (N, N) array_like
        Similarity array from SNF

    Returns
    -------
    W : (N, N) np.ndarray
        Stable-normalized similiarity array
    """

    # add `alpha` to the diagonal and symmetrize `W`
    rowSum = np.sum(W, 1) - np.diag(W)
    rowSum[rowSum == 0] = 1

    W = W / (2 * rowSum)
    np.fill_diagonal(W, 0.5)
    W = check_symmetric(W, raise_warning=False)

    return W