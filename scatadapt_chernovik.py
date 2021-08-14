import numpy as np
import scipy.stats as ss

def Ji(th,it,model=None, D=1):
    it = np.array([it]) if len(it.shape) == 1 else it
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
        return res
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
        return res2
    it = np.array([it]) if len(it.shape) == 1 else it
    if type(1)==int:
        th=np.array([th])
    if np.absolute(np.mean(it[:,1])-1)<0.00001 and np.var(it[:,1])<0.00001:
        mod="1PL"
    else:
        mod<-"other"
    data = dataGen(it.shape[0],model=mod)
    if (mod=="1PL"):
        res=np.empty((it.shape[0],1+len(th)))
        res[:]=np.nan
        for i in range(data.shape[0]):
            res[i, 0] = thetaEst(it, data[i], method = method, priorDist = priorDist,priorPar = priorPar, weight = weight,tuCo = tuCo, rang = rang)
        for j in range(len(th)):
            res[:, j+1]= LW(th[j], it)[:,1]
    else:
        res=np.empty((data.shape[0],1+len(th)))
        res[:]=np.nan
        for i in range(data.shape[0]):
            res[i, 0] = thetaEst(it, data[i], method = method, priorDist = priorDist,priorPar = priorPar, weight = weight,tuCo = tuCo, rang = rang)
            for j in range(len(th)):
                pi = Pi(th[j], it)['pi']
                qi = 1 - pi
                res[i, 1 + j] = np.prod(pi**data[i] * qi**(1 - data[i]))
    return res     



def Pi(th : float, it : np.array, model = None, D = 1):
    it = np.array([it]) if len(it.shape) == 1 else it 
    if model is None:
        a = it[:, 0]
        b = it[:, 1]
        c=  it[:,2]
        d=  it[:,3]
        e = np.exp(D * a * (th - b))
        pi = c + (d - c) * e/(1 + e)
        pi = np.where(pi == 0, 1e-10, pi)
        pi = np.where(pi == 1, 1 - 1e-10, pi)
        dpi = D * a * e /(1 + e)**2
        d2pi = D**2 * a**2 * e * (1 - e) /(1 + e)**3
        d3pi = D**3 * a**3 * e * (e**2 - 4 * e + 1)/(1 + e)**4
        res = {'pi' : pi, 'dpi' : dpi, 'd2pi' : d2pi, 'd3pi' : d3pi}
    return res
 
def integrate_cat(x, y):
    x1 = np.array(x[1:len(x)])
    x2 = np.array(x[0:len(x) - 1])
    hauteur = x1 - x2
    m = np.array([y[0:len(y) - 1], y[1:len(y)]])
    length = len(y)
    base = []
    i = 0
    for i in range(len(y) - 1):
        base.append((m[0, i] + m[1, i]) / 2)
        i += 1
    base = np.array(base)
    res = np.sum(hauteur * base)
    return res
  
def Ii(th, it, model = None, D = 1):
    it = np.array([it]) if len(it.shape) == 1 else it
    pr = Pi(th, it, model = model, D = D)
    P = pr['pi']
    dP = pr['dpi']
    d2P = pr['d2pi']
    d3P = pr['d3pi']
    if model is None:
        Q = 1 - P
        Ii = dP ** 2 / (P * Q)
        dIi = dP * (2 * P * Q * d2P - dP ** 2 * (Q - P)) / (P ** 2 * Q ** 2)
        d2Ii = (2 * P * Q * (d2P ** 2 + dP * d3P) - 2 * dP ** 2 * d2P * (Q - P)) / (P ** 2 * Q ** 2) - (3 * P ** 2 * Q * dP ** 2 * d2P - P * dP ** 4 * (2 * Q - P)) / (P ** 2 * Q ** 4)
    else:
        raise Exception('Не реализовано!')
    res = {'Ii': Ii, 'dIi': dIi, 'd2Ii': d2Ii}
    return res


def eapEst (it,x,model=None,D=1,priorDist="norm",priorPar=[0,1],lower=-4,upper=4,nqp=33):
    if model is None:
        def L(th,it,x):
            return np.prod(Pi(th,it,D=D)['pi']**x*(1-Pi(th,it,D=D)['pi'])**(1-x))
        def g(s):
            res=[]
            for i in range(len(s)):
                if priorDist=='norm':
                    res.append(s[i]*ss.norm.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='unif':
                    res.append(s[i]*ss.uniform.pdf(s[i],priorPar[1],priorPar[2])*L(s[i],it,x))
                elif priorDist=='Jeffreys':
                    res.append(s[i]*np.sqrt(sum(Ii(s[i],it,D=D)['Ii']))*L(s[i],it,x))
            return np.array(res)
        def h(s):
            res=[]
            for i in range(len(s)):
                if priorDist=='norm':
                    res.append(ss.norm.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='unif':
                    res.append(ss.uniform.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='Jeffreys':
                    res[i].append(np.sqrt(sum(Ii(s[i],it,D=D)['Ii']))*L(s[i],it,x))
            return np.array(res)
        X=np.linspace(start=lower,stop=upper,num=nqp)
        Y1=g(X)
        Y2=h(X)
    else:
        raise Exception('Не реализовано!')
    RES=integrate_cat(X,Y1)/integrate_cat(X,Y2)
    return RES

def eapSem (thEst,it,x,model=None,D=1,priorDist="norm",priorPar=[0,1],lower=-4,upper=4,nqp=33):
    if model is None:
        def L(th,it,x):
            return np.prod(Pi(th,it,D=D)['pi']**x*(1-Pi(th,it,D=D)['pi'])**(1-x))
        def g(s):
            res=[]
            for i in range(len(s)):
                if priorDist=='norm':
                    res.append((s[i]-thEst)**2*ss.norm.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='unif':
                    res.append((s[i]-thEst)**2*ss.uniform.pdf(s[i],priorPar[1],priorPar[2])*L(s[i],it,x))
                elif priorDist=='Jeffreys':
                    res.append((s[i]-thEst)**2*np.sqrt(sum(Ii(s[i],it,D=D)['Ii']))*L(s[i],it,x))
            return np.array(res)
        def h(s):
            res=[]
            for i in range(len(s)):
                if priorDist=='norm':
                    res.append(ss.norm.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='unif':
                    res.append(ss.uniform.pdf(s[i],priorPar[0],priorPar[1])*L(s[i],it,x))
                elif priorDist=='Jeffreys':
                    res[i].append(np.sqrt(sum(Ii(s[i],it,D=D)['Ii']))*L(s[i],it,x))
            return np.array(res)
        X=np.linspace(start=lower,stop=upper,num=nqp)
        Y1=g(X)
        Y2=h(X)
    else:
        raise Exception('Не реализовано!')
    RES=np.sqrt(integrate_cat(X,Y1)/integrate_cat(X,Y2))
    return RES

def nextItem (itemBank, model = None, theta = 0, out = [], x = [], 
    criterion = "MFI", method = "BM", priorDist = "norm", priorPar = [0,1], D = 1, rang = [-4, 4], parInt = [-4, 4, 33], infoType = "observed", 
    randomesque = 1, random_seed = None, rule = "length", thr = 20, 
    SETH = None, AP = 1, nAvailable = None, maxItems = 50, cbControl = None, 
    cbGroup = None):
    crit =criterion
    if cbControl is None:
        OUT = out
    else:
        raise Exception('Не реализовано!')
    if nAvailable is not None:
        raise Exception('Не реализовано!')
    if crit == "MFI":
        items = np.array([1 for i in range(itemBank.shape[0])] )
        for i in OUT:
            items[i] = 0
        info = Ii(theta, itemBank, model = model, D = D)['Ii']
        ranks=len(info)+1 - ss.rankdata(info).astype(int)
        nrIt = np.min(np.array([randomesque, sum(items)]))
        keepRank = np.sort(ranks[np.where(items == 1)])[:nrIt]
        keep =[]
        for i in range(len(keepRank)):
            p=np.where(ranks == keepRank[i])
            if items[p]==1:
                keep.append(p)
        if random_seed is not None:
            np.random.seed(random_seed)
        if len(keep)==1:
            select=keep[0]
        else:
            select=np.random.choice(keep, 1)
        res = {'item': select, 'par' : itemBank[select], 
            'info' : info[select], 'criterion': criterion, 'randomesque' : randomesque,'name':None}
    
    if crit == "bOpt":
        raise Exception('Не реализовано!')
    if crit == "MLWI" or crit == "MPWI":
        raise Exception('Не реализовано!')
    if crit == "KL" or crit == "KLP":
        raise Exception('Не реализовано!')
    if crit == "GDI" or crit == "GDIP":
        raise Exception('Не реализовано!')
    if crit == "MEI":
        raise Exception('Не реализовано!')
    if crit == "MEPV":
        raise Exception('Не реализовано!')
    if crit == "random":
        raise Exception('Не реализовано!')
    if crit == "progressive":
        raise Exception('Не реализовано!')
    if crit == "proportional":
        raise Exception('Не реализовано!')
    if crit == "thOpt":
        raise Exception('Не реализовано!')
    if cbControl is None:
        pass
    else:
        raise Exception('Не реализовано!')
    np.random.seed(None)
    return res
def startItems(itemBank, model = None, fixItems = None, nrItems = 1, theta = 0, D = 1, randomesque = 0, 
               seed = None, startSelect = "MFI", random_seed=None, nAvailable = None, cbControl = None, cbGroup = None, random_cb = None):
    if nAvailable is not None:
        raise Exception('Не реализовано!')
    if cbControl is not None:
        raise Exception('Не реализовано!')
    if nrItems > 0:
        if startSelect == "progressive" or startSelect == "proportional":
            raise Exception('Не реализовано!')
    if fixItems is not None:
        raise Exception('Не реализовано!')
    else:
        if nrItems > 0:
            if seed is not None: 
                np.random.seed(seed)
                if cbControl is None:
                    if nAvailable is None: 
                        items = np.random.choice(np.array([i for i in range(itemBank.shape[0])]), nrItems)
                    else:
                        items = np.random.choice(np.where(nAvailable == 1),   nrItems)
                else:
                    np.random.seed(None)
                    raise Exception('Не реализовано!')
                par = itemBank[items]
                startSelect = None
                thStart =  None
                np.random.seed(None)
            else:
                if type(theta) == int:
                    theta=[theta]
                thStart = theta
                if startSelect == "bOpt":
                    raise Exception('Не реализовано!')
                if startSelect == "thOpt":
                    raise Exception('Не реализовано!')
                if startSelect == "MFI":
                    items = np.array([])
                    nr_items = itemBank.shape[0]
                    selected =np.tile(0, nr_items)
                    if random_seed is not None:
                        np.random.seed(random_seed)
                    for i in range(len(thStart)):
                        item_info = Ii(thStart[i], itemBank, model = model, D = D)['Ii']
                        if nAvailable is not None: 
                            pos_adm = (1 - selected) * np.array(nAvailable)
                        else:
                            pos_adm = 1 - selected
                        pr = np.sort(item_info[np.where(pos_adm == 1)])
                        keep = pr[randomesque]
                        k = np.where(item_info >= keep)
                        prov=np.array([])
                        for y in k[0]:
                            if pos_adm[y]==1:
                                prov=np.append(prov,y)
                        if len(prov) == 1:
                            items=np.append(items, int(prov[0]))
                        else:
                            items=np.append(items, int(np.random.choice(prov, 1)[0]))
                        selected[int(items[0])] = 1
                items=items.astype(int)
                par = itemBank[items]
            res = {'items' :items, 'par' : par, 'thStart' : thStart, 
                'startSelect' : startSelect,'names':None}
        else:
            res = {'items': None, 'par': None, 'thStart': None, 
            'startSelect': None,'names': None}
    np.random.seed(None)
    return(res)

def checkStopRule(th, se, N, it = None, model = None, D=1, stop):
    res = False
    res_rule = None
    for i in range(len(stop['rule'])):
        if stop['rule'][i] == 'length':
            ind = 1
        if stop['rule'][i] == 'precision':
            ind = 2
        if stop['rule'][i] == 'classification':
            ind = 3
        if stop['rule'][i] == 'minInfo':
            ind = 4
        if ind == 1:
            if N >= stop['thr'][i]:
                res = True
                res_rule = np.array([res_rule, stop['rule'][i])
        if ind == 2:
            if se <= stop['thr'][i]:
                res = True
                res_rule = np.array(res_rule, stop['rule'][i])
        if ind == 3:
            if th - ss.norm.ppf(1-stop['alpha']/2) * se >= stop['thr'][i] or th + ss.norm.ppf(1-stop['alpha']/2) * se <= stop['thr'][i] :
                res = True
                res_rule = np.array(res_rule, stop['rule'][i])
        if ind == 4:
            info = Ii(th, it, model = model, D = D)['Ii']
            if max(info) <= stop['thr'][i]:
                res = True
                res_rule = np.array(res_rule, stop['rule'][i])
    return {'desicion': res, 'rule': res_rule}
