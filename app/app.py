import streamlit as st
import pandas as pd
import numpy as np


st.title('Fecundability Predictions')


chart_data = pd.DataFrame(np.random.randn(20, 3), columns=["a", "b", "c"])
st.bar_chart(chart_data)
