import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = []
y = []
z = []

with open("partData") as file:
    lines = file.readlines()
    for line in lines[1:]:
        elems = line.split(' ')
        x.append(float(elems[0]))
        y.append(float(elems[1]))
        z.append(float(elems[2]))

#plt.xlim(-10000, 10000)
#plt.ylim(-10000, 10000)
#plt.zlim(-10000, 10000)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label='parametric curve')
#ax.set_xlim(-10000, 10000)
#ax.set_ylim(-10000, 10000)
#ax.set_zlim(-10000, 10000)
plt.show()

#x_old_arr = []
#y_old_arr = []
#x_new_arr = []
#y_new_arr = []

#with open("partData") as file:
#    lines = file.readlines()
#    for line in lines[0:]:
#        elems = line.split(' ')
#        x_old_arr.append(float(elems[0]))
#        y_old_arr.append(float(elems[1]))
#        x_new_arr.append(float(elems[2]))
#        y_new_arr.append(float(elems[3]))
#plt.plot(x_old_arr, y_old_arr)
#print("Ploting...");
#plt.show()
