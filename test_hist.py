import numpy as np
import matplotlib.pyplot as plt
import os
# Fixing random state for reproducibility
np.random.seed(19680801)

mu, sigma = 100, 15
x = mu + sigma * np.random.randn(10000)

# the histogram of the data
#n, bins, patches = plt.hist(x, 50, density=True, outlinecolor='Blue',facecolor='darkturquoise', alpha=0.75)
n, bins, patches = plt.hist(x, 50, density=True, ec='Blue',color='darkturquoise', alpha=0.75)

n, bins, patches = plt.hist(x, 50, density=True, edgecolor='Blue',color='darkturquoise', alpha=0.75,histtype= "stepfilled", linewidth=1.5)
plt.style.use('dark_background')
plt.show()


plt.hist(x,50)
style=plt.style.available
for st in style:
    plt.figure(st)
    plt.xlabel('Smarts')
    plt.ylabel('Probability')
    plt.title('Histogram of IQ')
    plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    plt.xlim(40, 160)
#plt.ylim(0, 0.03)
    plt.grid(True)

    print(st)
    plt.style.use(st)
    plt.hist(x,50)
    plt.show()

#plt.show()
