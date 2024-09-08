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


chart_data = aneuploid_prob()
st.bar_chart(chart_data, x="Maternal Age", y_label="Probability")
