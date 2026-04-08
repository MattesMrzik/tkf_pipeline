import numpy as np

def gap_concentration(df):
    df['gap_concentration'] = np.where(
        df['gap_col%'] > 0,
        df['gap%'] / df['gap_col%'],
        0
    )
    return df
