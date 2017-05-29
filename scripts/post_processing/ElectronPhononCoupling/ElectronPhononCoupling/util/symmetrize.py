import numpy as np

def symmetrize_array(arr, indices_groups, axis=0):
    """
    Symmetrize an array by taking the average over each subset of indices.

    For example, given a list of indices_groups [[a1,a2,a3], [b1,b2,b3]],
    the arrays (arr[a1], arr[a2], arr[a3]) will be averaged,
    as well as the arrays (arr[b1], arr[b2], arr[b3]).

    Arguments
    ---------

    arr: np.array
        The array to be symmetrized.
    indices_groups: 2D list.
        A list of groups of indices
    axis:
        The dimension on which the indicies are taken.

    Returns
    -------

    sym_arr: np.array
        The symmetrized array, of the same shape as arr.
    """

    sym_arr = np.copy(arr)

    indices = np.arange(arr.size).reshape(arr.shape)

    for group in indices_groups:

        arr_group = arr.take(group, axis=axis)
        av_arr = np.average(arr_group, axis=axis) 

        for i in group:

            ind = indices.take(i, axis=axis)
            sym_arr.put(ind, av_arr)

    return sym_arr



