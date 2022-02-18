# NeuralNetwork
BasicAnalysis and number of motifs 
#Code_Motifs_Matteo


from numpy import genfromtxt
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random as ran
import pandas as pd
import itertools
    



def transform_To_Undirected_Network (connectivity_matrix):
    """
    transform the directed network connectivity matrix to an undirected one
    removing also the self-loop.

    Parameters
    ----------
    connectivity_matrix : 2D array
        connectivity matrix of the directed network

    Returns
    -------
    new_matrix : 2D array
        connectivity matrix of the undirected equivalent network

    """
    new_matrix = connectivity_matrix
    a,b = new_matrix.shape
    new_matrix[np.identity(a)==1] = 0
    temp = np.transpose(new_matrix)
    new_matrix = (new_matrix + temp)>=1 
    new_matrix = new_matrix.astype('uint8')
    return new_matrix

def NumberEdges(connectivity_matrix):
    return np.sum(connectivity_matrix>0)




def density_Network(connectivity_matrix, network_type = "directed"):
    """
    compute 

    Parameters
    ----------
    connectivity_matrix : TYPE
        DESCRIPTION.
    network_type : TYPE, optional
        DESCRIPTION. The default is "directed".

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    a, b = connectivity_matrix.shape
    matrix = connectivity_matrix >= 1
    matrix = matrix.astype('uint8')
    if network_type == "directed" :
        density = np.sum(matrix)/(a*a)
        return density
    elif network_type == "undirected":
        
        density = (np.sum(matrix)/2)/((a*(a-1))/2)
        return density
    else :
        print ("you entered the wrong type of network !")
        return 0


def inOut_degree (connectivity_matrix):
    a, b = connectivity_matrix.shape
    degree_in = np.zeros((a))
    degree_out = np.zeros((a))
    matrix = connectivity_matrix >= 1
    
    matrix = matrix.astype('uint8')
    matrix[np.identity(a) == 1] = 0 #on ne compte pas les self loop
    # print(matrix)
    for line in range(a):
        for column in range(a):
            degree_out [column] += matrix[column,line]
            degree_in[line] += matrix[column,line]
    return degree_in, degree_out


        
def SCC_detector (connectivity_matrix):
    """
    Return List of Nodes for each SCC

    Parameters
    ----------
    connectivity_matrix : 2D array
        connectivity matrix of the directed network

    Returns
    -------
    SCC : list of array
        list of nodes for each SCC
"""
    
    
    # first we need to know link that goes in 2 ways
    a, b = connectivity_matrix.shape
    unweighted = np.zeros((a,b))
    unweighted[connectivity_matrix != 0] = 1
    reciprocalEdge = np.transpose(unweighted) + unweighted == 2
    reciprocalEdge = reciprocalEdge.astype('uint8')
    reciprocalEdge [np.identity(a) == 1] = 0
    
    # Then we do a lists of nodes connecting to each other (connected to i node element) it function as sub SCC






def twoD_motifNumber(connectivity_matrix):
    notweighted_connectivity_matrix=connectivity_matrix
    notweighted_connectivity_matrix[connectivity_matrix>0]= 1
    T = notweighted_connectivity_matrix + np.transpose(notweighted_connectivity_matrix)
    
    twod1=0
    twod2=0
    N=np.array(np.shape(connectivity_matrix))[0]

    
    for i in range(1,N):
        for j in range (0,i):
           
            if T[i,j]==1:
                twod1 += 1
            if T[i,j]==2:
                twod2 += 1
               
    return twod1, twod2

def threeD_motifNumber(connectivity_matrix):
    #count all the possible motifs of 3 nodes, the first three definition 
    #are the same we gave for the twoD_motifNumber
    M=connectivity_matrix
    M[connectivity_matrix>0]= 1
    T = M + np.transpose(M)
    N=np.array(np.shape(connectivity_matrix))[0]
    
    #below xe define for each different kind of 3 nodes motifs a variable that will count how many 
    #to explain we will use this notation: -> directional link, -><- double link, i,j,k are our three nodes.
    threed1=0 
    #  i -> j, i -> k and simmetric ones
    threed2=0
    # i -> j, j -> k
    threed3=0
    # 
    threed4=0
    threed5=0
    threed6=0
    threed7=0
    threed8=0
    threed9=0
    threed10=0
    threed11=0
    threed12=0
    threed13=0
    
    for i in range(1,N):
        for j in range (i):
            for k in range (j):
                a=T[i,j];
                b=T[j,k];
                c=T[k,i];
                
                #in this way we are computing the total number of links in our motif. We are considering 3 nodes motif with at least 2 links. So we start counting the ones with 2 links:
                if (a+b+c)==2:
                    
                    # Here we found the difference between the first 
                    if (M[i,j]==1 and M[i,k]==1) or (M[j,k]==1 and M[j,i]==1) or (M[k,i]==1 and M[k,j]==1):
                        threed1 +=1;
                        
                    if (M[i,j]==1 and M[j,k]==1) or (M[k,i]==1 and M[i,j]==1) or (M[i,k]==1 and M[k,j]==1) or (M[k,j]==1 and M[j,i]==1) or (M[i,k]==1 and M[j,i]==1) or (M[k,i]==1 and M[j,k]==1):
                        threed3 +=1;
                        
                    if (M[i,j]==1 and M[k,j]==1) or (M[i,k]==1 and M[j,k]==1) or (M[j,i]==1 and M[k,i]==1):
                        threed2 +=1;
                
                # all the 3 nodes motifs with 3 links    
                if (a+b+c)==3:
                    
                    if a==2:
                        # 3 nodes motifs with a double link between i and j
                        if M[i,k]==1 or M[j,k]==1:
                            threed4 +=1;
                        
                        if M[k,i]==1 or M[k,j]==1:
                            threed5 +=1;
                            
                    if b==2:
                        #double link between j and k
                        if M[j,i]==1 or M[k,i]==1:
                            threed4 +=1
                            
                        if M[i,j]==1 or M[i,k]==1:
                            threed5 +=1;
                            
                    if c==2 :
                        #double link between k and i 
                        if M[k,j]==1 or M[i,j]==1:
                            threed4 += 1;
                            
                        if M[j,k]==1 or M[j,i]==1:
                            threed5 +=1;
                    
                    # three single link
                    if a==1 and b==1 and c==1:
                        
                        if (M[i,j]==1 and M[j,k]==1 and M[k,i]==1) or (M[i,k]==1 and M[k,j]==1 and M[j,i]==1):
                            threed6 +=1
                            # the links creating a loop
                        if (M[i,j]==1 and M[k,j]==1 and M[k,i]==1) or (M[i,j]==1 and M[k,j]==1 and M[i,k]==1) or (M[i,j]==1 and M[j,k]==1 and M[i,k]==1) or (M[j,i]==1 and M[j,k]==1 and M[i,k]==1) or (M[j,i]==1 and M[j,k]==1 and M[k,i]==1) or (M[j,i]==1 and M[k,j]==1 and M[k,i]==1):
                            threed7 +=1
                            # 2 links in one direction, the other in the opposite one    
                if (a+b+c)==4:
                    if a==0 or b==0 or c==0:
                        threed10 +=1;
                        
                    if a==2:
                        if M[j,k]==1 and M[i,k]==1:
                            threed9 +=1
                            
                        if M[k,i]==1 and M[k,j]==1:
                            threed11 +=1
                            
                        if M[j,k]==1 and M[k,i]==1: 
                            threed8 +=1
                            
                    if b==2:
                        if M[j,i]==1 and M[k,i]==1:
                            threed9 +=1
                            
                        if M[i,j]==1 and M[i,k]==1:
                            threed11 +=1
                            
                        if M[k,i]==1 and M[i,j]==1:
                            threed8 +=1
                            
                    if c==2:
                        if M[i,j]==1 and M[k,j]==1:
                            threed9 +=1
                            
                        if M[j,i]==1 and M[j,k]==1:
                            threed11 +=1
                            
                        if M[i,j]==1 and M[j,k]==1:
                            threed8 +=1
                            
                if a+b+c==5:
                    threed12 +=1
                    
                if a+b+c==6:
                    threed13 +=1
                        
                
    return [threed1, threed2, threed3, threed4, threed5, threed6, threed7, threed8, threed9, threed10, threed11, threed12, threed13]


#Null models comparison

def ErdosRenyi (Nodes, p1):
    connectivityNull=np.zeros((Nodes, Nodes))
    n=0
    for i in range(Nodes):
        for j in range(Nodes):
            if ran.random() <= p1:
                if connectivityNull[i,j]==0:
                    connectivityNull[i,j]=1
    return connectivityNull
        
def edges_list(N):
    
    n=np.array(np.shape(N))[0]
    
    x=[]
    for i in range(n):
        for j in range(n):
            if N[i,j]!=0:
                x.append((i,j))
    
    return x, np.matrix(x)

def nullModelDegreeD(arr,n):
    
    edges_in=arr[:,0]
    edges_out=arr[:,1]
    print(arr[:,0])
    nn=len(arr[:,0])
    print(nn)
    D=np.zeros([n,n])
    m=np.arange(nn)
    q=np.arange(nn)

    i=0
    print(m)
    while i < nn:
        
        b=np.random.choice(m)
        np.delete(m,np.where(m==b))
        a=np.random.choice(q)
        np.delete(q,np.where(q==a))
        
        if D[int(edges_in[a]),int(edges_out[b])] == 0:
            D[edges_in[a],edges_out[b]]=1
            i+=1
    return D


def Network_from_matrix(N,edges,method=None):
    if method==None:
        G=nx.empty_graph(N,create_using=nx.Graph)
    else:
        G.empty_graph(N,create_using=method)
    G.add_edges_from(edges)
    return G

def K_core(G,K):
    NODES=np.zeros([len(list(G.nodes())),len(list(G.nodes()))])
    nodes=list(G.nodes())
    core_di=[]
    removed=[]
    NODES[0]=nodes
    for i in range(len(nodes)):
        if G.degree[i]<K:
            nodes.remove(i)
            removed.append(i)
            
    
    Gk=G.subgraph(nodes)
    NODES[1,0:len(nodes)]=nodes
    k=0
    while np.array_equal(NODES[k,:],NODES[k+1,:],equal_nan=False)==False:
        for i in list(Gk.nodes()):
            if Gk.degree[i]<K:
                nodes.remove(i)
                removed.append(i)
        
        NODES[k+2,0:len(nodes)]=nodes
        k+=1
        Gk=Gk.subgraph(nodes)
    core_di.append(len(removed))
    n=len(list(Gk.nodes()))
    return Gk, n, nodes, core_di
        
    
        
    
        
        
            
    
def IN_Strenght_Distribution(connect):
    in_strenght=[]
    for i in range(np.shape(connect)[0]):
        in_strenght.append(sum(connect[i,:]))
    return in_strenght

def OUT_Strenght_Distribution(connect):
    out_strenght=[]
    for i in range(np.shape(connect)[0]):
        out_strenght.append(sum(connect[:,i]))
    return out_strenght

def strenght_distribution(connect):
    strenght=[]
    for i in range(np.shape(connect)[0]):
        strenght.append(sum(connect[i,:])+sum(connect[:,i]))
    return strenght

def weight_distribution(connect):
    weight=[]
    for i in range(np.shape(connect)[0]):
        for j in range(np.shape(connect)[0]):
            if connect[i,j]!=0:
                weight.append(connect[i,j])
    return weight



  
        
                    
#%% process

brutData = genfromtxt(r"C:\Users\matte\Downloads\connectivity matrix table1 project.csv", delimiter=',')

brutData2 =  genfromtxt(r"C:\Users\matte\Downloads\connectivity matrix table1 project.csv", delimiter=',', dtype=None, encoding='utf-8')


N=387

connect = brutData[1:N+1, 1:N+1]

undirected = transform_To_Undirected_Network(connect)

d1 = density_Network(connect,network_type = "directed")

d2 = density_Network(connect, network_type = "undirected")

a, b = inOut_degree(connect)

a = np.array(a)

b = np.array(b)

edges=NumberEdges(connect)

strenght=strenght_distribution(connect)
plt.hist(strenght,100)
plt.title("Strenght distribution")
plt.xlabel("strenght")
plt.ylabel("number of nodes")
plt.show()
plt.close()

weight=weight_distribution(connect)
plt.hist(weight,70)
plt.xlabel("weight")
plt.ylabel("number of links")
plt.title("Weight distribution")
plt.show()
plt.close()

# undirected_edge(40,mydata)
NumberOfEdges = NumberEdges(connect) 



x,xx=edges_list(connect)

G=Network_from_matrix(N,x)

#the colormap, to plot the network
colormap=[]
MBON=[]
OAN=[]
DAN=[]
APL=[]
MBIN=[]
PN=[]
KC=[]

for i in range(1,N+1):
    if "MBON" in str(brutData2[i,0]):
        colormap.append("blue")
        MBON.append(i-1)
    if "OAN" in str(brutData2[i,0]):
        colormap.append("pink")
        OAN.append(i-1)
    if "DAN" in str(brutData2[i,0]):
        colormap.append("orange")
        DAN.append(i-1)
    if "APL" in str(brutData2[i,0]):
        colormap.append("cyan")
        APL.append(i-1)
    if "MBIN" in str(brutData2[i,0]):
        colormap.append("yellow")
        MBIN.append(i-1)
    if "PN" in str(brutData2[i,0]):
        colormap.append("green")
        PN.append(i-1)
    if "KC" in str(brutData2[i,0]):
        colormap.append("red")    
        KC.append(i-1)
        
colorplot=["blue","pink","orange","cyan","yellow","green","red"]
y=np.zeros([7,68])

k=0
nodes=list(G.nodes())
Nodes=[]
pos=nx.spring_layout(G,k=0.8  )
nx.draw(G,pos,node_size=10 , node_color=colormap, node_shape='o', alpha=1.0)
plt.title("Our network")
plt.show()
core_dist=[]
Gk, Nk, nodes, core_d = K_core(G, k)
core_dist.append(N-core_d[0])
Nodes.append(nodes)

while nodes!=[]:
    k=k+1
    Gk, Nk, nodes, core_d = K_core(G, k)
    core_dist.append(N-core_d[0])
    Nodes.append(nodes)
    y[0,k-1]=(len(set(nodes).intersection(set(MBON)))/len(MBON))
    y[1,k-1]=(len(set(nodes).intersection(set(OAN)))/len(OAN))
    y[2,k-1]=(len(set(nodes).intersection(set(DAN)))/len(DAN))
    y[3,k-1]=(len(set(nodes).intersection(set(APL)))/len(APL))
    y[4,k-1]=(len(set(nodes).intersection(set(MBIN)))/len(MBIN))
    y[5,k-1]=(len(set(nodes).intersection(set(PN)))/len(PN))
    y[6,k-1]=(len(set(nodes).intersection(set(KC)))/len(KC))
    nx.draw_networkx_nodes(G,pos,node_size=10 , node_color=colormap, node_shape='o', alpha=1.0)
    nx.draw_networkx_edges(G, pos, Gk.edges())
    plt.title("K_core, K="+ str(k))
    plt.show()
plt.close()

Gk,Nk,nodes,a = K_core(G, 67)
Nodes.append(nodes)
nx.draw_networkx_nodes(G,pos,node_size=10 , node_color=colormap, node_shape='o', alpha=1.0)
nx.draw_networkx_edges(G, pos, Gk.edges())
plt.title("K_core, K=67")
plt.show()
plt.close()
print(xx[:,0])

plt.plot(core_dist[0:len(core_dist)-1], "ro")
plt.ylabel("K core size")
plt.xlabel("K")
plt.show()
plt.close()

for i in range(7):
    plt.plot((range(len(y[i,:])-1)),y[i,0:len(y[i,:])-1],color=colorplot[i])
plt.title("Core filling factor")
plt.ylabel("fraction of nodes in the K core by type")
plt.xlabel("K")
plt.show()
plt.close()

degree_in, degree_out = inOut_degree(connect)

y2=twoD_motifNumber(connect)
plt.hist(y2)
plt.show()
plt.close()


y3=threeD_motifNumber(connect)
q=np.arange(1,len(y3)+1)
plt.plot(q,y3)
plt.title("Neural network connection patterns")
plt.ylabel("number of motifs")
plt.xlabel("type of motif")
plt.show()
plt.close()

#theoretical result
n=len(list(itertools.combinations(list(G.nodes),3)))
p=d2
y1=[n*3*p**2*(1-p)**4,n*3*p**2*(1-p)**4,n*6*p**2*(1-p)**4,n*6*p**3*(1-p)**3,n*6*p**3*(1-p)**3,n*2*p**3*(1-p)**3,n*6*p**3*(1-p)**3,n*6*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*6*p**5*(1-p),n*p**6]
plt.plot(range(1,14),y1, "-r")
plt.title("Expected result")
plt.ylabel("number of motifs")
plt.xlabel("type of motif")
plt.show()
plt.close()
p=d2
n=len(list(itertools.combinations(range(N),3)))
y1=[n*3*p**2*(1-p)**4,n*3*p**2*(1-p)**4,n*6*p**2*(1-p)**4,n*6*p**3*(1-p)**3,n*6*p**3*(1-p)**3,n*2*p**3*(1-p)**3,n*6*p**3*(1-p)**3,n*6*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*3*p**4*(1-p)**2,n*6*p**5*(1-p),n*p**6]
plt.plot(range(1,14),y1, "-r")
plt.title("Erdos Renyi graph N=387, p=d=" + str(d2))
plt.ylabel("number of motifs")
plt.xlabel("type of motif")
m=np.zeros([13,10])

for i in range(10):
    connectivityNullD=ErdosRenyi(N,d2) 
    y = threeD_motifNumber(connectivityNullD)
    for j in range(13):
        m[j,i]=y[j]
    q=np.arange(1,len(y)+1)
    plt.plot(q,y)
plt.legend(["expected result"])
plt.show()
plt.close()
    
for i in range(13):
    plt.hist((m[i][:]))
    plt.title("m="+str(i))
    plt.show()
    plt.close()

m=np.zeros([13,10])

for i in range(10):
    connectivityNullD=nullModelDegreeD(xx,N)
    y = threeD_motifNumber(connectivityNullD)
    for j in range(13):
        m[j,i]=y[j]
    q=np.arange(1,len(y)+1)
    plt.plot(q,y)
plt.title("Random model preserving in and out degree of the nodes")
plt.ylabel("number of motifs")
plt.xlabel("type of motif")
plt.show()
plt.close()
    
for i in range(13):
    plt.hist((m[i][:]))
    plt.title("m="+str(i))
    plt.show()
    plt.close()

        
