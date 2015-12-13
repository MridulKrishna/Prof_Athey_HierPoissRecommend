#Script to plot the likelihood versus the iteration
f1 = open('output.txt','r')
lines = f1.readlines()
wide_core = []
wide_L2size = []
narrow_core = []
narrow_L2size = []
for line in lines:
	index = line.find("Likelihood = ")
	wide_core = wide_core + [int(line1[index + 7 : index + 8])]
	index = line1.find("Size = ")
	wide_L2size = wide_L2size + [int(float(line1[index + 7:]))]
for line2 in lines2:
	index = line2.find("Core = ")
	narrow_core = narrow_core + [int(line2[index + 7 : index + 9])]
	index = line2.find("Size = ")
	narrow_L2size = narrow_L2size + [int(float(line2[index + 7:]))]
print wide_core
print wide_L2size
print narrow_core
print narrow_L2size
import matplotlib.pyplot as plt
plt.figure(0)
plt.bar(wide_core,wide_L2size)
plt.xlabel('Cores')
plt.ylabel('Max L2 Cache Sizes (MB)')
plt.savefig('Max_L2cacheSizes_wide.png')
plt.figure(1)
plt.bar(narrow_core,narrow_L2size)
plt.xlabel('Cores')
plt.ylabel('Max L2 Cache Sizes (MB)')
plt.savefig('Max_L2cacheSizes_narrow.png')

