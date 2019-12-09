import numpy as np
import random
from numpy import linalg as la
from matplotlib import pyplot as plt 

n_real = 100


def recur_factorial(n):  
   if n == 1:  
       return n  
   else:  

       return n * recur_factorial(n-1) 
n_dataP = 20
q = 0
VV = np.zeros(n_dataP)
R_avg = np.zeros(n_real)
RF_avg = np.zeros(n_dataP)
for V_o in np.linspace(0.5, 10.0, n_dataP):
    
    VV[q] = V_o

    for o in range(n_real):
        #if f == 0:
            #n_spins = 10
       # else:
        n_spins = 10
        gamma = 0.
        J = 2.    
        V = V_o * np.random.rand(n_spins - 1)    
        
        Alpha = 0.
        count = 0
        
        n_ham = int( recur_factorial(n_spins)/ ( recur_factorial(n_spins/2)**2) ) #number of basis states that are used after conditioning with 0 magnetisation
    
        #print("nham" + str(n_ham))
        Ham = np.zeros((2**n_spins, 2**n_spins))
        s_a = np.zeros(n_ham).astype(int) #array of numbers whose binary representation indicates spin arrangement in chain
    
        index = 0
        for i in range(2**n_spins):
            count = 0
            for ch in bin(i):
                if ch == '1':
                    count+= 1 
            if count == n_spins/2: ## checking if half of the sites are occupied
                s_a[index] = int(i)
                index += 1
    
                
    
        def findstate(sx): ##searching for state sx in array s_a
            
            bmin = 0
            bmax = n_ham
            while(True):
                b = int(bmin + (bmax - bmin)/2)
    
                if (sx < s_a[b]):
                    bmax = b - 1
                elif (sx > s_a[b]):
                    bmin = b + 1
                else:
                    return b
    
        # def findstate2(sx): 
        #         return list(s_a).index(sx)
    
        for a in range(n_ham): 
            for i in range(n_spins - 1):
              
              j = i + 1
    
              if (format((s_a[a]),'0'+str(n_spins)+'b')[i] == format((s_a[a]),'0'+str(n_spins)+'b')[j]): # if adjacent spins are same or not 
    
                  Ham[s_a[a]][s_a[a]] = Ham[s_a[a]][s_a[a]] + 1/4 * V[i]
              else:
                
                  Ham[s_a[a]][s_a[a]] = Ham[s_a[a]][s_a[a]] -1/4 * V[i]
                  
                  sx = (s_a[a] ^ (1<<(n_spins - 1 - i) ) )^ ( 1<<(n_spins - 1 - j) ) 
                  
                  b = findstate(sx) 
            
                  Ham[s_a[a]][s_a[b]] = Ham[s_a[a]][s_a[b]]+ 1/2 * J 
    
            for i in range(n_spins):
              W = -gamma * (i) + Alpha * (((i)/ (n_spins-1 ))**2)
              if (format(( s_a[a]), '0' + str(n_spins) + 'b' )[i] == '1'): #site at ith position has spin down or up
                  Ham[s_a[a]][s_a[a]] = Ham[s_a[a]][s_a[a]] - W /2
              elif(format(( s_a[a]), '0' + str(n_spins) + 'b' )[i] == '0'):
                  Ham[s_a[a]][s_a[a]] = Ham[s_a[a]][s_a[a]] + W /2
    
        #print(Ham)
    
        Eval, Evec = la.eig(Ham)
        Einterest = []
        for x in range(len(Eval)):
          if Eval[x] != 0.:
            Einterest.append(Eval[x]) #adding non zero elements of list of eigenvalues to a new array called Einterest
          
        Eval, Evec = la.eig(Ham)
        #print( "Eval {}".format(Eval))
        
        #print( len(Einterest) == n_ham)
        Einterest= np.array(Einterest)
        Einterest.sort()
        #print(Einterest)
        r = np.zeros(((len(Einterest) - 2)))
        sum = 0
        for i in range(len(Einterest) - 2):
          d1 = Einterest[i+1] - Einterest[i]
          d2 = Einterest[i+2] - Einterest[i+1]
          if (d1 == 0 or d2 == 0):
            r[i] = 0
          r[i] = min( ((d1/d2), (d2/d1)))
          
        #print( "r.mean {}".format(r.mean()))
        R_avg[o] = r.mean()
        #print("V {}".format(V))
        #print("J {}".format(J))
    RF_avg[q] = R_avg.mean()
    q = q+1 
    

#print( R_avg)
if( n_spins == 8):
    col = "or"
else:
    col = "xb"

plt.title("Title ({} , {} )".format(n_spins, n_real)) 
plt.xlabel("V") 
plt.ylabel("R_AVG") 
plt.plot(VV, RF_avg, 'xr')
plt.show()


#plt.plot(VV, R_avg, col)
#plt.show() 



  
            
    
# print( "i" +str( i) )
      # print(bin(s_a[a][0]))
      # sy = s_a[a][0] ^ (1<<(n_spins - 1 - i))
      # print(bin(sy))
      # b = findstate(sy)
      # if ( format((s_a[a][0]),'0'+str(n_spins)+'b')[i]) =='0':
      #   Ham[s_a[a][0]][s_a[a][0]] = Ham[s_a[a][0]][s_a[a][0]] + 0.5 * W
      # else:
      #   Ham[s_a[a][0]][s_a[a][0]] = Ham[s_a[a][0]][s_a[a][0]] - 0.5 * W

        


