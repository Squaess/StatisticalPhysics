import numpy as np
import matplotlib.pyplot as plt
import csv
from sklearn.linear_model import LinearRegression

graph = []

with open('data.csv') as f:
    csv_reader = csv.reader(f, delimiter=",")
    for row in csv_reader:
        graph.append((float(row[0]), float(row[1])))

X = np.log(np.array([x[0] for x in graph])).reshape(-1,1)
y = np.log(np.array([x[1] for x in graph])).reshape(-1,1)

reg = LinearRegression().fit(X, y)
print(reg.coef_)
domain = np.arange(0, 10**5, 0.1).reshape(-1,1)
prediction = reg.predict(domain)
plt.plot(domain, prediction)
plt.xlabel("log(Numbers of steps)")
plt.ylabel("log(variance)")
plt.show()