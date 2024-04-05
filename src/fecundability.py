import numpy as np


class FecundabilityEstimator:
    """Create an estimator of"""

    def __init__(self):
        pass

    def p_aneuploid(self):
        """Return the probability of a specific aneuploidy type.

        Entries correspond to baseline probability of:
        * meiotic origin
        * tripolar-mitotic
        * low mosaic
        * high mosaic
        * euploid
        """
        return np.array([0.4, 0.05, 0.1, 0.05, 0.4])

    def p_failed_implantation(self):
        """The probability of failed implantation conditional on aneuploidy type.

        Entries correspond to the probability of a failed implantation conditional on:
        * meiotic origin aneuploidy
        * tripolar-mitotic
        * low mosaic
        * high mosaic
        * euploid
        """
        return np.array([0.8, 1.0, 0.4, 0.5, 0.05])

    def p_epl_implantation(self, eta=1e-2, eps=1e-1):
        """Probability of early pregnancy loss conditional on implantation occurring.

        The  probability of an early-pregnancy loss conditional on implantation succeeding

        The output vector corresponds to the probability:
        * meiotic origin aneuploidy
        * tripolar-mitotic
        * low mosaic
        * high mosaic
        * euploid

        Params:
            * eta (`float`): This is the "escape probability" of a meiotic or tripolar mitotic aneuploidy.
            * eps (`float`): This is the probability of euploid embryo EPL
        """
        assert eta < 1.0
        assert eps < 1.0
        return np.array([1.0 - eta, 1.0 - eta, 0.2, 0.5, eps])
