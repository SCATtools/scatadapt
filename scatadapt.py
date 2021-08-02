import numpy as np
def Ji(th,it,model=None, D=1):
  pr=Pi(th,it,model=model,D=D)
  P = pr['pi']
  dP = pr['dpi']
  d2P = pr['d2pi']
  d3P = pr['d3pi']
  if model is None:
    Q = 1-P
    ji = dP * d2P/(P * Q)
    dji = (P * Q * (d2P**2 + dP * d3P) - dP**2 * d2P *(Q - P))/(P**2 * Q**2)
    res = {'ji': ji, 'dji' : dji}
  else:
    raise Exception('Не реализовано!')
  return res
def fullDist(th, it, method = "BM", priorDist="norm",priorPar=np.array([0,1]), weight = "Huber", tuCo   = 1, range  = np.array([-4 ,4]), parInt = np.array([-4, 4, 33])):
  def dataGen(n, model="1PL"):
    if model=="1PL":
      res=np.zeros((n+1,n))
      for i in range(1,n+1):
        res[i , 0:i] = 1
    else:
      res = numpy.empty((2**n,n))
      for i in range(n):
        k=0
        t=0
        for y in range(2**n):
          if k==2**(n-i-1):
            k=0
            t=(t+1)%2
          res[y,i]=t
          k+=1
    return(res)
  def LW(th,it,D=1):
    P=Pi(th,it,D=D)['pi']
    Q=1-P
    res=np.empty((it.shape[0]+1,it.shape[0]))
    res[:]=np.nan
    res[0,0]=Q[0]
    res[1,0]=P[0]
    for i in range(1,it.shape[0]):
      for j in range(i+2):
        if j==0:
          res[j,i]=res[j,i-1]*Q[i]
        else:
          if j==i:
            res[j,i]<-res[j-1,i-1]*P[i]
              else:
                res[j,i]<-res[j-1,i-1]*P[i]+res[j,i-1]*Q[i]
    res2=np.concatenate((np.reshape(np.array(range(res.shape[1]+1)),(res.shape[1]+1,1)),np.reshape(res[:,-1],(res.shape[1]+1,1))),axis=1)
    return(res2)
             
def Pi():
  
