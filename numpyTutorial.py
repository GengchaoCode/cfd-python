import numpy as np

# The use of numpy built-in Numpy functions can provide an increase in execution speed many-times over
u = np.array([0, 1, 2, 3, 4, 5])

for i in range(1, len(u)):
    print(u[i]-u[i-1])

# calculate using the numpy way via slicing
print(u[1:]-u[:-1])
