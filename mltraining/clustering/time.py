import numpy as np

from matplotlib import pyplot as plt
import seaborn as sns

fname = "time.dat"

with open(fname, "r") as file:

    times = []

    for line in file:
        time = line.split()[0]

        m, s = time.split(":")

        times.append(60 * float(m) + float(s))

times = np.array(times) / 60

print(len(times), np.mean(times), np.std(times))

sns.boxplot(y=times)
plt.ylabel("Time (min)")
plt.show()