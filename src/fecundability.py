import numpy as np
import matplotlib.pyplot as plt
import warnings


class FecundabilityEstimator:
    """Estimators for the conditional probabilities involved in pregnancy losses."""

    def __init__(self):
        pass

    def p_aneuploid(self, X=None):
        """Return the probability of a specific aneuploidy type.

        Entries correspond to baseline probability of:
        ```
             | Meiotic | Euploid
        None | X1      |  X2
        Early| X3      |  X4
        Late | X5      |  X6
        ```
        """
        if X is None:
            return np.array(
                [
                    [294 / 909, 206 / 909],
                    [260 / 909 * 3 / 4, (88 + 61) / 909 * 3 / 4],
                    [260 / 909 * 1 / 4, (88 + 61) / 909 * 1 / 4],
                ]
            )
        else:
            assert np.all((X > 0.0) & (X < 1.0))
            assert np.isclose(X.sum(), 1.0)
            return X

    def scale_meiotic_prob(self, A, p=0.2):
        """Scale the probability of meiotic aneuploidy."""
        assert (p > 0) & (p < 1.0)
        assert A.shape[0] == 3
        assert A.shape[1] == 2
        X = A.copy()
        # We have to scale everything to the original margianl expected probability of aneuploidy ...
        p_m = X[:, 0].sum()
        X[:, 0] = X[:, 0] * (p / p_m)
        X[:, 1] = X[:, 1] * ((1 - p) / (1 - p_m))
        assert np.isclose(X.sum(), 1.0)
        return X

    def p_failed_implantation(self, X=None):
        """The probability of failed implantation conditional on aneuploidy type.

        Entries correspond to the probability of a failed implantation conditional on:
        ```
             | Meiotic | Euploid
        None | X1      |  X2
        Early| X3      |  X4
        Late | X5      |  X6
        ```
        """
        if X is None:
            return np.array([[0.36, 0.16], [0.55, 0.57], [0.36, 0.16]])
        else:
            assert np.all((X > 0.0) & (X < 1.0))
            return X

    def p_epl_implantation(self, eta=1e-2, eps=1e-1):
        """Probability of early pregnancy loss conditional on implantation occurring.

        The  probability of an early-pregnancy loss conditional on implantation succeeding

        Entries correspond to the probability of EPL (biochemical pregnancy) conditional on:
        ```
             | Meiotic | Euploid
        None | X1      |  X2
        Early| X3      |  X4
        Late | X5      |  X6
        ```
        Params:
            * eta (`float`): This is the "escape probability" of a meiotic or tripolar mitotic aneuploidy.
            * eps (`float`): This is the probability of euploid embryo EPL
        """
        assert eta < 1.0
        assert eps < 1.0
        return np.array([[1.0 - eta, eps], [1.0 - eta, 1.0 - eta], [1.0 - eta, eps]])

    def p_miscarriage(self, eta_m=1e-2, eps_m=1e-1):
        """Probability of miscarriage conditional on not EPL & implantation not failing

        Entries correspond to the probability of Miscarriage conditional on:
        ```
             | Meiotic | Euploid
        None | X1      |  X2
        Early| X3      |  X4
        Late | X5      |  X6
        ```
        Params:
            * eta_m (`float`): This is the "escape probability" of a meiotic or tripolar mitotic aneuploidy.
            * eps_m (`float`): This is the probability of euploid embryo EPL

        """
        return np.array(
            [[1.0 - eta_m, eps_m], [1.0 - eta_m, 1.0 - eta_m], [1.0 - eta_m, eps_m]]
        )

    def p_implant_marginal(self, A=None, A_implant=None, gamma_implant=5e-2):
        """Marginal probability of an implantation error."""
        assert A.shape == (3, 2)
        p_implant = gamma_implant + np.sum(
            self.p_failed_implantation(X=A_implant) * self.p_aneuploid(X=A)
        )
        if p_implant > 1.0:
            warnings.warn("Marginal probability of implantation error is > 1!")
        return p_implant

    def p_epl_marginal(self, A=None, gamma_epl=5e-2, gamma_implant=5e-2, **kwargs):
        """Marginal probability of an early pregnancy loss."""
        assert A.shape == (3, 2)
        assert (gamma_epl > 0) and (gamma_epl < 1.0)
        p_epl = gamma_epl + np.sum(
            self.p_epl_implantation(**kwargs)
            * (1.0 - self.p_failed_implantation() - gamma_implant)
            * self.p_aneuploid(X=A)
        )
        if p_epl > 1.0:
            warnings.warn("Marginal probability of EPL is > 1!")
        return p_epl

    def p_miscarriage_marginal(
        self,
        A=None,
        A_implant=None,
        gamma_mis=5e-2,
        gamma_epl=5e-2,
        gamma_implant=5e-2, **kwargs,
    ):
        """Marginal probability of miscarriage."""
        assert A.shape == (3, 2)
        assert (gamma_mis > 0) and (gamma_mis < 1.0)
        assert (gamma_epl > 0) and (gamma_epl < 1.0)
        assert (gamma_implant > 0) and (gamma_implant < 1.0)
        p_mis = gamma_mis + np.sum(
            self.p_miscarriage(**kwargs)
            * (1.0 - self.p_epl_implantation(**kwargs) - gamma_epl)
            * (1.0 - self.p_failed_implantation(X=A_implant) - gamma_implant)
            * self.p_aneuploid(X=A)
        )
        if p_mis > 1.0:
            warnings.warn("Marginal probability of miscarriage is > 1!")
        return p_mis


class FecundabilityPlotting:
    def __init__(self):
        pass

    def plot_macklon_pyramid(self, ax, xs=np.array([0.3, 0.3, 0.1, 0.3]), **kwargs):
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
        txts = []
        for i, x in enumerate(xs):
            txts.append(ax.text(0.5, pos[i], f"{x*100:0.1f}%", **kwargs))
        return ax, txts
