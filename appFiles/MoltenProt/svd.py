import numpy as np


def apply_svd(X):
    """
    Perform Singular Value Decomposition (SVD) on the input data matrix X.

    Parameters
    ----------
        X : numpy array of shape (n_wavelengths, n_measurements)
            The input data matrix to decompose.

    Returns
    -------
        explained_variance : numpy array
            The cumulative explained variance for each component.
        basis_spectra     : numpy array
            The left singular vectors (U matrix) representing the basis spectra.
        coefficients      : numpy array
            The coefficients associated with each basis spectrum.
    """

    U, S, _ = np.linalg.svd(X)

    # Calculate the total variance or correlation
    total_variance = np.sum(S ** 2)
    cumulative_variance = np.cumsum(S ** 2)

    # The matrix V contains the variation of each component against the temperature / measurement dimension

    a_is = []

    for i in range(U.shape[1]):
        def coefficients_bi(column):
            return U[:, i].dot(column)

        a_i = np.apply_along_axis(coefficients_bi, axis=0, arr=X)

        a_is.append(a_i)

    coefficients = np.array(a_is)

    # Basis spectra
    basis_spectra = U

    # Cumulated explained variance of the components
    explained_variance = cumulative_variance / total_variance * 100

    return explained_variance, basis_spectra, coefficients