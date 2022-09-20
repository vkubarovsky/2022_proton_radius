import numpy as np
import matplotlib.pyplot as plt
  


# let's make a line model
def line(x, a, b):
    return a + x * b

a_true = 1.0
b_true = 2.0

a_fit=1.0
b_fit=2.0

# let's make some data
x = np.linspace(0, 1, 10)
xx=np.array([0,1.,1.,0.5])
# precomputed random numbers from standard normal distribution
z = np.array([-0.49783783, -0.33041722, -1.71800806,  1.60229399,
                 1.36682387, -1.15424221, -0.91425267, -0.03395604,
                 -1.27611719, -0.7004073 ])

sigma_y = 0.1 * np.ones_like(x)
y = line(x, a_true, b_true) + sigma_y * z

#plt.errorbar(x, y, sigma_y, fmt="o")
plt.xlim(-0.1, 1.1);

#plt.errorbar(x, y, sigma_y, fmt="o")
plt.plot(xx, line(xx, a_fit, b_fit))
plt.xlim(-0.1, 1.1);

plt.show()

