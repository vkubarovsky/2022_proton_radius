from iminuit import Minuit

def cost_function(x, y, z):
    return (x - 2) ** 2 + (y - 3) ** 2 + (z - 4) ** 2

m = Minuit(cost_function, x=0, y=0, z=0)

m.migrad()  # run optimiser
m.hesse()   # run covariance estimator

print(m.values)  # x: 2, y: 3, z: 4
print(m.errors)  # x: 1, y: 1, z: 1
