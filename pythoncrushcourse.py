# <-- comments in Python are denoted by the pound sign, like this one

## Libraries
import numpy as np                      # we import the array library
from matplotlib import pyplot as plt    # import the plotting library (pyplot is a module in the matplotlib library)

# Each function name is written following the library name, with a dot in between
myarray = np.linspace(1,5,11)
print(myarray)

## Variables
# Python does not require explicitly declared variable type like C and other languages
a = 5                                   # a is an integer 5
b = 'five'                              # b is a string of the word 'five
c = 5.0                                 # c is a floating point 5

print(type(a))
print(type(b))
print(type(c))

## Whitespace in Python
# Python uses indents and whitespace to group statements together
'''
for (i=0, i<5, i++){
    printf('Hi! \n');
}
'''
for i in range(5):
    print('Hi! \n')

# If you have nested for-loops, there is a further indent for the inner loop
for i in range(3):                      # start from 0, upto but not include 3 (3 loops in total)
    for j in range(3):
        print(i, j)

    print('This statement is within the i-loop, but not the j-loop.')

## Slicing Arrays
# In Python, you can look at portions of arrays in the same way as in Matlab, with a few extra tricks thrown in.
myvals = np.array([1, 2, 3, 4, 5])
print(myvals)

# Python uses zero-based index, so the first and the last element in the array myvals is:
print(myvals[0], myvals[4])

# The slice in Python is inclusive in the front and exclusive in the end
print(myvals[0:3])

## Assigning Array Variables
# The equal sign in Python just create an alias, instead of a copy (assignment by reference) --> save memory
a = np.linspace(1,5,5)
b = a
a[2] = 17
print(a)
print(b)
print(id(a))
print(id(b))

# Make a true copy of a
c = a.copy()
a[2] = 3
print(a)
print(c)
print(id(a))
print(id(c))