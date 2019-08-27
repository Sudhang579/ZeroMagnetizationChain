
import numpy as np
from numpy import linalg as la

def recur_factorial(n):  
   if n == 1:  
       return n  
   else:  
       return n*recur_factorial(n-1) 

gamma = 0.
J = 2.
V = 1.
Alpha = 1.
count = 0
n_spins = 12

n_ham = int( recur_factorial(n_spins)/ ( recur_factorial(n_spins/2)**2) ) #number of basis states that are used after conditioning with 0 magnetisation

print("n_ham" + str(n_ham))
Ham = np.zeros((2**n_spins, 2**n_spins)) 
s_a = np.zeros(n_ham).astype(int) #matrix of numbers 

index = 0
for i in range(2**n_spins):
    count = 0
    for ch in bin(i):
        if ch == '1':
            count+= 1
    if count == n_spins/2:
        s_a[index] = int(i)
        index += 1

        

def findstate(sx):
    
    bmin = 1
    bmax = n_ham 
    while(True):
        b = int(bmin + (bmax - bmin)/2)

        if (sx < s_a[b]):
            bmax = b - 1
        elif (sx > s_a[b]):
            bmin = b + 1
        else:
            return b

def findstate2(sx):
        return list(s_a).index(sx)

for a in range(n_ham):
    for i in range(n_spins - 1):
      
      j = i + 1

      if (format((s_a[a]),'0'+str(n_spins)+'b')[i] == format((s_a[a]),'0'+str(n_spins)+'b')[j]):

          Ham[a][a] = Ham[a][a] + 1/4 * V
      else:
         
          Ham[a][a] = Ham[a][a] -1/4 * V
          
          sx = (s_a[a] ^ (1<<(n_spins - 1 - i) ) )^ ( 1<<(n_spins - 1 - j) ) 
          
          b = findstate2(sx) 
     
          Ham[a][b] = Ham[a][b]+ 1/2 * J ###

    for i in range(n_spins):
      W = -gamma * (i) + Alpha * (((i)/ (n_spins-1))**2)
      if (format(( s_a[a]), '0' + str(n_spins) + 'b' )[i] == '1'):
          Ham[a][a] = Ham[a][a] - W /2
      elif(format(( s_a[a]), '0' + str(n_spins) + 'b' )[i] == '0'):
          Ham[a][a] = Ham[a][a] + W /2

print(Ham)

Eval, Evec = la.eig(Ham)
Einterest = [0]
for x in range(len(Eval)):
  if Eval[x] != 0.:
    Einterest.append(Eval[x])
print( len(Einterest) == n_ham)
Einterest= np.array(Einterest)
Einterest.sort()
print(Einterest)
r = np.zeros(((len(Einterest) - 2)))
sum = 0
for i in range(len(Einterest) - 2):
  d1 = Einterest[i+1] - Einterest[i]
  d2 = Einterest[i+2] - Einterest[i+1]
  if (d1 == 0 or d2 == 0):
    r[i] = 0
  r[i] = min( ((d1/d2), (d2/d1)))
  
print( r.mean())






  
            
    
# print( "i" +str( i) )
      # print(bin(s_a[a][0]))
      # sy = s_a[a][0] ^ (1<<(n_spins - 1 - i))
      # print(bin(sy))
      # b = findstate(sy)
      # if ( format((s_a[a][0]),'0'+str(n_spins)+'b')[i]) =='0':
      #   Ham[s_a[a][0]][s_a[a][0]] = Ham[s_a[a][0]][s_a[a][0]] + 0.5 * W
      # else:
      #   Ham[s_a[a][0]][s_a[a][0]] = Ham[s_a[a][0]][s_a[a][0]] - 0.5 * W

        


