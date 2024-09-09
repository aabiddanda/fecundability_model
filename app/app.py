import streamlit as st
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

from libfecund import FecundabilityEstimator

st.title("Fecundability Predictions")

gruhn2019_df = pd.read_csv("../data/gruhn_etal2019.fig1E.csv")

pmeiotic_age_day3 = interp1d(
    gruhn2019_df[gruhn2019_df.sample_type_indicator == 0].maternal_age.values,
    gruhn2019_df[gruhn2019_df.sample_type_indicator == 0].fitted_values.values,
    fill_value="extrapolate",
)

pmeiotic_age_day5 = interp1d(
    gruhn2019_df[gruhn2019_df.sample_type_indicator == 1].maternal_age.values,
    gruhn2019_df[gruhn2019_df.sample_type_indicator == 1].fitted_values.values,
    fill_value="extrapolate",
)

f_est = FecundabilityEstimator()


def aneuploid_prob(maternal_age=np.arange(20, 52), day3=True):
    if day3:
        ps = pmeiotic_age_day3(maternal_age)
    else:
        ps = pmeiotic_age_day5(maternal_age)
    A = f_est.p_aneuploid()
    As = [f_est.scale_meiotic_prob(A, p=a) for a in ps]
    data = np.array(
        [
            [As[i][0, 0] for i in range(ps.size)],
            [As[i][1, 0] for i in range(ps.size)],
            [As[i][2, 0] for i in range(ps.size)],
            [As[i][0, 1] for i in range(ps.size)],
            [As[i][1, 1] for i in range(ps.size)],
            [As[i][2, 1] for i in range(ps.size)],
        ]
    )
    data_df = pd.DataFrame(
        {
            "Maternal Age": maternal_age,
            "Meiotic-only": data[0, :],
            "Meiotic & Early-Mitotic": data[1, :],
            "Meiotic & Late-Mitotic": data[2, :],
            "Euploid-only": data[3, :],
            "Euploid & Early-Mitotic": data[4, :],
            "Euploid & Late-Mitotic": data[5, :],
        }
    )
    return data_df


def pregnancy_outcome(maternal_age=np.arange(20, 52), day3=True):
    if day3:
        ps = pmeiotic_age_day3(maternal_age)
    else:
        ps = pmeiotic_age_day5(maternal_age)
    A = f_est.p_aneuploid()
    As = [f_est.scale_meiotic_prob(A, p=a) for a in ps]
    outcome_probs = []
    for A in As:
        p_mis, p_epl, p_implant = (
            f_est.p_miscarriage_marginal(A=A),
            f_est.p_epl_marginal(A=A),
            f_est.p_implant_marginal(A=A),
        )
        # Should we have a re-scaling here to ensure sum to 1?
        if (p_mis + p_epl + p_implant) > 1.0:
            x = np.array([0, p_mis, p_epl, p_implant])
        else:
            x = np.array([1.0 - (p_mis + p_epl + p_implant), p_mis, p_epl, p_implant])
        x /= np.sum(x)
        outcome_probs.append(x)

    outcome_probs = np.vstack(outcome_probs)
    data_df = pd.DataFrame(
        {
            "Maternal Age": maternal_age,
            "Implantation Error": outcome_probs[:, 3],
            "Early-Pregnancy Loss": outcome_probs[:, 2],
            "Miscarriage": outcome_probs[:, 1],
            "Live Birth": outcome_probs[:, 0],
        }
    )
    return data_df


min_age = st.slider("Minimum Maternal Age", 15, 55, 18)
max_age = st.slider("Maximum Maternal Age", 15, 55, 52)
st.subheader("Aneuploidy Occurrence")
aneuploidy_data = aneuploid_prob(maternal_age=np.linspace(min_age, max_age, 100))
st.area_chart(aneuploidy_data, x="Maternal Age", y_label="Probability", stack=True)

st.subheader("Pregnancy Losses")
outcome_data = pregnancy_outcome(maternal_age=np.linspace(min_age, max_age, 100))
st.area_chart(outcome_data, x="Maternal Age", y_label="Probability", stack=True)
