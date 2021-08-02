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
             
def Pi():
  
