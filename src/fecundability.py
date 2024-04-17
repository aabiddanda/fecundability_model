import numpy as np
import matplotlib.pyplot as plt


class FecundabilityEstimator:
    """Estimators for the conditional probabilities involved in pregnancy losses."""

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


class FecundabilityPlotting:
    def __init__():
        pass

    def plot_macklon_pyramid(ax, xs=np.array([0.3, 0.3, 0.1, 0.3]), **kwargs):
        """Plot a version of the Macklon et al pyramid."""
        assert xs.size == 4
        assert np.isclose(xs.sum(), 1.0)
        ax.plot([0, 0.5], [0, 1], color="black")
        ax.plot([0.5, 1], [1, 0], color="black")
        ax.plot([0, 1], [0, 0], color="black")
        # Divide pyramid into four tiers
        ax.plot([0.25 / 2, 1.0 - 0.25 / 2], [0.25, 0.25], lw=2, color="black")
        ax.plot([0.5 / 2, 1.0 - 0.5 / 2], [0.5, 0.5], lw=2, color="black")
        ax.plot([0.75 / 2, 1.0 - 0.75 / 2], [0.75, 0.75], lw=2, color="black")
        ax.axhline(0.5, lw=0.5, linestyle="--", color="black")
        pos = [0.25 / 2, (0.25 + 0.5) / 2, (0.5 + 0.75) / 2, (0.75 + 0.9) / 2]
        for i, x in enumerate(xs):
            ax.text(0.5, pos[i], f"{x*100}%", **kwargs)
        return ax
