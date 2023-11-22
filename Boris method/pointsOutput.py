import matplotlib.pyplot as plt

x_old_arr = []
y_old_arr = []
x_new_arr = []
y_new_arr = []

with open("partData") as file:
    lines = file.readlines()
    for line in lines[0:]:
        elems = line.split(' ')
        x_old_arr.append(float(elems[0]))
        y_old_arr.append(float(elems[1]))
        x_new_arr.append(float(elems[2]))
        y_new_arr.append(float(elems[3]))
plt.plot(x_old_arr, y_old_arr)
print("Ploting...");
plt.show()
