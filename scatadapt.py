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
    prov=dP*d2P/P
    prov1=(P*d2P**2+P*dP*d3P-dP**2*d2P)/P**2
    res=as.numeric(rowSums(prov,na.rm=TRUE))
    resd=as.numeric(rowSums(prov1,na.rm=TRUE))
    res={Ji:res,dJi:resd}
  return res
def Pi():
  
