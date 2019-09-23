import numpy as np

BGO = np.loadtxt('Compton_110deg_15min_ch000.txt', skiprows=5, usecols=(0,1))
NaI = np.loadtxt('Compton_110deg_15min_ch001.txt', skiprows=5, usecols=(0,1))
ns = 10**(-9) #1 nanosecond

lim_BGOt = np.array([])   #Initiating new array to be able to append in for loop. (Optimization reasons)
lim_NaIt = np.array([])
lim_BGOc = np.array([])
lim_NaIc = np.array([])

#Simplified version,but slow. Made to gain an insight of the more complex loop.
#for j in Nai:
#    for i in BGO:
#       d= j[0] - i[0]
#        if( d < 5 and d > -5):
#             lim_BGOt = np.append(lim_BGOt, i[0])
#             lim_NaIt = np.append(lim_NaIt, j[0])
#             lim_BGOc = np.append(lim_BGOc, i[1])
#             lim_NaIc = np.append(lim_NaIc, j[1])

#Optimized version of the simpler
#The longest list must be in the outer loop.
for i in sorted(NaI[:,0]):        #Sorting NaI and BGO
    for j in sorted(BGO[:, 0]):
            d = j - i  # Taking the time difference of their entrances
            if (d*ns > 100*ns): #Checking wether the coincidence difference is larger than a set timelimit
                # Since NaI and BGO is sorted we make use of that we know that next element will be larger
                break
            elif (d*ns < -100*ns):
                # If the difference is negative with larger difference than the timelimit keep searching through the list
                continue
            #elif( 700 > np.abs(BGO[np.where(BGO == j)[0][0]][1] + NaI[np.where(NaI == i)[0][0]][1]) > 600):
            elif(True):

                #The elif tests for energy conservation between the scattered photon and recoil electron
                #Set true if energy conservation is neglected (Remember to calibrate before if taking EC into account).

                lim_BGOt = np.append(lim_BGOt, j)  #Appending entrance to array
                lim_NaIt = np.append(lim_NaIt, i)  #Nai time list
                lim_BGOc = np.append(lim_BGOc, BGO[np.where(BGO == j)[0][0]][1])   #Appending to BGO channel list
                lim_NaIc = np.append(lim_NaIc, NaI[np.where(NaI == i)[0][0]][1])   #Appending to BGO time list

#Combining the time and channel lists into a single list of 2 coloums and N number of rows.
lim_BGO = np.transpose((lim_BGOt,lim_BGOc))
lim_NaI = np.transpose((lim_NaIt,lim_NaIc))

#Saving the data into .txt files, 1 couloumn for time, 1 for channelnumber
np.savetxt('filtered_BGO_only_Time.txt', lim_BGO, delimiter=" ", fmt='%s')
np.savetxt('filtered_NaI_only_Time.txt', lim_NaI, delimiter=" ", fmt='%s')

