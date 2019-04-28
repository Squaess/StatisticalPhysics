import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv('flips.csv')
plt.scatter(df['k'], df['m'])
plt.show()