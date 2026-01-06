import pandas as pd

df1 = pd.read_csv("WH_gen_sample_summary.csv")
df2 = pd.read_csv("ZH_gen_sample_summary.csv")
df3 = pd.read_csv("ggH_gen_sample_summary.csv")

merged = pd.concat([df1, df2, df3], ignore_index=True)
merged.to_csv("HtoBBGG_gen_sample_summary.csv", index=False)

