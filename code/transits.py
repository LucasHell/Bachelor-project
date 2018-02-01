import matplotlib.pyplot as plt

timesFile = []
valueArray = []
transitTime1Float = []
epoch1Float = []
count = 0


timesFile = open('Times')
valueArray = timesFile.readlines()
planet = [k.split(' ')[0] for k in valueArray]
epoch1 = [k.split(' ')[1] for k in valueArray]
transitTime1 = [k.split(' ')[2] for k in valueArray]


for k in range(len(valueArray)):
	if planet[k] == '0':
		epoch1Float.append(float(epoch1[k]))
		transitTime1Float.append(float(transitTime1[k]))
		count += 1


plt.plot(epoch1Float, transitTime1Float, label='Transit Time')
plt.xlabel('Epoch')
plt.ylabel('Transit time')
plt.title('Transit Timing variations')
plt.legend()
plt.show()

