import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('Mass_Radius.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=';')
    for row in plots:
        y.append(float(row[0]))
        x.append(float(row[1]))

plt.plot(x,y, label='Mass of exoplanet')
plt.xlabel('Planet Radius $R_E$')
plt.ylabel('Planet Mass $M_E$')
plt.title('Mass-radius relation')
plt.legend()
plt.show()
