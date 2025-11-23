import pandas as pd

df1 = pd.read_csv("WH_gen_sample_summary.csv")
df2 = pd.read_csv("ZH_gen_sample_summary.csv")

merged = pd.concat([df1, df2], ignore_index=True)
merged.to_csv("VH_gen_sample_summary.csv", index=False)

