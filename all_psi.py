# Mu Jinglin,牟景林
#ref1,https://github.com/bryanjp3/OrbitalViewer/blob/master/wave.py
#ref2,https://github.com/rafaariza/Orbitales/blob/master/orbitales.py
"""
需要预先安装的模块
python -m pip install numpy -i https://pypi.doubanio.com/simple/
python -m pip install mayavi -i https://pypi.doubanio.com/simple/
python -m pip install matplotlib -i https://pypi.doubanio.com/simple/
"""
# Importing libraries in alphabetical order
from datetime import datetime  # Library for working with dates and times
import math  # Library for mathematical functions
from math import factorial
import os  # Library for interacting with the operating system
os.environ['ETS_TOOLKIT'] = 'qt4'
#添加代码,不显示warnings
import warnings
warnings.filterwarnings('ignore')
#代码结束

import matplotlib  # Library for creating visualizations
import matplotlib.pyplot as plt  # Sublibrary of matplotlib for creating plots
from mpl_toolkits.mplot3d import Axes3D
import numpy as np  # Library for working with arrays

from mayavi import mlab  # Library for creating 3D visualizations
from mayavi.core.api import PipelineBase  # Base class for all Mayavi pipelines
from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel  # Classes for creating Mayavi scenes

from traits.api import HasTraits, Instance, Range, on_trait_change  # API for defining classes with traits
from traitsui.api import Group, Item, View  # API for creating graphical user interfaces

#添加代码,正常显示中文
from pylab import mpl 
mpl.rcParams['font.sans-serif'] = ['SimHei'] 
plt.rcParams['axes.unicode_minus']=False
#代码结束

pi=np.pi
a = 5.291772108e-11

def input_nlm():
    #input n
    while True :
        try:
            n = int(input("1st quantum number n: "))
            if n <= 0:
                print('n应该大于0')
                continue
            else:
                break
        except:
            print('n应该是整数')
    #input l
    while True :
        try:
            l = int(input("2nd quantum number l: "))
            if l < 0 or l>=n:
                print('l应该大于等于0且小于n')
                continue
            else:
                break
        except:
            print('l应该是整数')
    if l < 0 or l>=n:
        print('l>0 and l<n')

    #input m
    while True :
        try:
            m = int(input("3rd quantum number m: "))
            if (m > l) or (m < -l):
                print('m应该属于[-l,l]')
                continue
            else:
                break
        except:
            print('l应该是整数')
    if m>l or m<(-l):
        print('Quantum number m belongs to [-l,l]')
    return n,l,m

def input_nl():
    #input n
    while True :
        try:
            n = int(input("1st quantum number n: "))
            if n <= 0:
                print('n应该大于0')
                continue
            else:
                break
        except:
            print('n应该是整数')
    #input l
    while True :
        try:
            l = int(input("2nd quantum number l: "))
            if l < 0 or l>=n:
                print('l应该大于等于0且小于n')
                continue
            else:
                break
        except:
            print('l应该是整数')
    if l < 0 or l>=n:
        print('l>0 and l<n')

    return n,l

def input_lm():
    #input l
    while True :
        try:
            l = int(input("quantum number l: "))
            if l < 0:
                print('l应该大于等于0')
                continue
            else:
                break
        except:
            print('l应该是整数')

    #input m
    while True :
        try:
            m = int(input("quantum number m: "))
            if (m > l) or (m < -l):
                print('m应该属于[-l,l]')
                continue
            else:
                break
        except:
            print('l应该是整数')
    if m>l or m<(-l):
        print('Quantum number m belongs to [-l,l]')
		
    return l,m

#计算legendre_polynomial
def legendre_polynomial(l,m,x):
    pmm = 1.0
    if m > 0:
        sign = 1.0 if m % 2 == 0 else -1.0
        pmm = sign*pow(factorial(2*m-1)*(1.0-x*x),((m/2)))

    if l == m:
        return pmm

    pmm1 = x*(2*m+1)*pmm
    if l == m+1:
        return pmm1

    for n in range(m+2,l+1):
        pmn = (x*(2*n-1)*pmm1-(n+m-1)*pmm)/(n-m)
        pmm = pmm1
        pmm1 = pmn

    return pmm1

#把r,Theta,Phi转为x,y,z
def sph2cart(r,Theta,Phi):

    x = r*np.sin(Theta)*np.cos(Phi)
    y = r*np.sin(Theta)*np.sin(Phi)
    z = r*np.cos(Theta)
    
    return x, y, z

#把x,y,z转为r,Theta,Phi
def cart2sph(x,y,z):
    phi = np.arctan2(y,x)
    theta = np.arctan2(np.sqrt(x**2 + y**2),z)
    r = np.sqrt(x**2 + y**2 + z**2)
    return r, theta, phi
    
#计算associatedLaguerre多项式
def associatedLaguerre(l, n, x):
    coeffs = []
    for k1 in range(n-l-1+1):
        c = (-1)**k1*factorial(n+l)**2/factorial(n-l-k1-1)/factorial(2*l+1+k1)/factorial(k1)
        coeffs.append(c)
        
    out = [0.0 for z in range(len(x))] # output
    for m in range(len(coeffs)):
        #print m, " ", coeffs[m], "jgfk"
        out = out + coeffs[m]*(x**m)
    return out

#计算波函数数值
def calc_psi(r,Theta,Phi,n, l, m):

    A = np.sqrt(((2*l+1)*factorial(l-abs(m)))/(4*pi*factorial(l+abs(m))))

    rho = 2.0*r/n/a
    
    if m>0:
            
        scalars = np.sqrt(2)*np.cos(m*Phi)*legendre_polynomial(l,m,np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
            
    elif m < 0:
            
        scalars = np.sqrt(2)*np.sin(abs(m)*Phi)*legendre_polynomial(l,abs(m),np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
            
    else:
            
        scalars = legendre_polynomial(l,0,np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
    
    #角度部分归一化常数
    A = np.sqrt(((2*l+1)*factorial(l-abs(m)))/(4*pi*factorial(l+abs(m))))
    #径向部分归一化常数
    B = np.sqrt((2.0/n)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3))
    #波函数归一化
    scalars = scalars * A * B

    return scalars

#计算波函数数值，并检查是否为0
def calc_psi_prof(r,Theta,Phi,n, l, m):

    
    rho = 2.0*r/n/a
    
    if m>0:
            
        scalars = np.sqrt(2)*np.cos(m*Phi)*legendre_polynomial(l,m,np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
            
    elif m < 0:
            
        scalars = np.sqrt(2)*np.sin(abs(m)*Phi)*legendre_polynomial(l,abs(m),np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
            
    else:
            
        scalars = legendre_polynomial(l,0,np.cos(Theta))*np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
    
    #角度部分归一化常数
    A = np.sqrt(((2*l+1)*factorial(l-abs(m)))/(4*pi*factorial(l+abs(m))))
    #径向部分归一化常数
    B = np.sqrt((2.0/n)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3))
    #波函数归一化
    scalars = scalars * A * B

    #波函数在某个平面是否取零,不能通过角度的方式屏蔽，原因在于可能出现零的零次方
    if np.max(abs(scalars)) < 0.01: scalars = scalars * 0

    return scalars

#计算球谐函数
def calc_Y(Theta,Phi,l, m):


    #print('cos(m*Phi)')
    #print(np.cos(m*Phi))
    #print('cos(Theta)')
    #print(np.cos(Theta))
    #print('sin(abs(m)*Phi)')
    #print(np.sin(abs(m)*Phi))
    if m>0:
        #print("m>0")
        scalars = np.sqrt(2)*np.cos(m*Phi)*legendre_polynomial(l,m,np.cos(Theta))

    elif m < 0:
        #print("m<0")
        scalars = np.sqrt(2)*np.sin(abs(m)*Phi)*legendre_polynomial(l,abs(m),np.cos(Theta))

    else:
        #print("m=0")
        scalars = legendre_polynomial(l,0,np.cos(Theta))

    #角度部分归一化常数
    A = np.sqrt(((2*l+1)*factorial(l-abs(m)))/(4*pi*factorial(l+abs(m))))
    #波函数归一化
    scalars = scalars * A 
    #波函数是否取零,不能通过角度的方式屏蔽，原因在于可能出现零的零次方
    if np.max(abs(scalars)) < 0.01: scalars = scalars * 0

    return scalars




def Psi_points():

    n,l,m = input_nlm()
    
    x, y, z = np.ogrid[-8*a*n:8*a*n:368j, -8*a*n:8*a*n:368j, -8*a*n:8*a*n:368j]
    ngrid=368
    r,Theta,Phi = cart2sph(x,y,z)
    
    scalars = calc_psi(r,Theta,Phi,n, l, m)
    #print('scalars')
    #print(scalars)
    #PlotNodes
    x1,y1,z1 = sph2cart(r,Theta,Phi) #再次转换，保持xyz和psi有合适的维度
    mlab.contour3d(x1, y1, z1, scalars, contours=[0], opacity=0.2, transparent=True, color=(0.5,0.5,0.0))
    sca_aver = np.average(abs(scalars))
    sca_aver = 500*sca_aver**2
    

    
    
    rand1 = np.random.rand(ngrid*ngrid*ngrid)-0.5
    rand2 = rand1.reshape([ngrid,ngrid,ngrid])
    x = x + rand2*a
    
    rand1 = np.random.rand(ngrid*ngrid*ngrid)-0.5
    rand2 = rand1.reshape([ngrid,ngrid,ngrid])
    y = y + rand2*a
    
    rand1 = np.random.rand(ngrid*ngrid*ngrid)-0.5
    rand2 = rand1.reshape([ngrid,ngrid,ngrid])
    z = z + rand2*a
    
    # 获取并显示第一个时间点
    time1 = datetime.now()
    time1_str = time1.strftime("%H:%M:%S")
    print(f"开始时间：{time1_str}")
	
    ###原来的耗时代码
    #arr = np.zeros((ngrid,ngrid,ngrid), dtype=bool)
	#
    #for i in range(len(scalars)):
    #    for j in range(len(scalars[i])):
    #            for k in range(len(scalars[i][j])):
    #                if (scalars[i][j][k])**2 < np.random.rand()*sca_aver:
    #                    arr[i][j][k]=False
    #                else:
    #                    arr[i][j][k]=True

    ###新代码
    rand_sca = np.random.rand(ngrid,ngrid,ngrid)
    arr = np.square(scalars) > sca_aver*rand_sca
    
	# 获取并显示第二个时间点
    time2 = datetime.now()
    time2_str = time2.strftime("%H:%M:%S")
    print(f"结束时间：{time2_str}")
    
    # 计算并显示时间差
    delta = time2 - time1
    print(f"运行时长：{delta}")
	
    scalars = scalars[arr]
    x = x[arr]
    y = y[arr]
    z = z[arr]
    scalars=scalars/abs(scalars)
    
    mlab.points3d(x, y, z, scalars, mode='point')#, colormap="blue-red"
    #mlab.colorbar()
    
    mlab.title('probability density cloud of Psi')
    mlab.show()




def Psi_superpos():


    class MyModel(HasTraits):
        posx = Range(0,10,10)
        posy = Range(0,10,0)
        posz = Range(0,10,0)
        coef1 = Range(-50,50,10)
        coef2 = Range(-50,50,10)
    	#将来可以添加只点击这个变化的按钮。
        #Draw = Button('Draw') 
    
    
        #场景模型实例
        scene = Instance(MlabSceneModel,()) #后面加上()是将他实例化了
        #管线实例
        plot = Instance(PipelineBase)
    
        def __init__(self,**traits):
            HasTraits.__init__(self,**traits)
    
            x2 = x-self.posx * a * n2**1.65
            y2 = y-self.posy * a * n2**1.65
            z2 = z-self.posz * a * n2**1.65
            r,Theta,Phi = cart2sph(x2,y2,z2)
    		
            scalars2 = calc_psi(r,Theta,Phi,n2,l2,m2)
    		
            scalars=scalars1*self.coef1 + scalars2*self.coef2
    
            maxpsi=np.abs(np.max(scalars))
            minpsi=np.abs(np.min(scalars))
            print("maxpsi",np.max(scalars))
            print("minpsi",np.min(scalars))
            limpsi=np.maximum(maxpsi,minpsi)
    
            if self.plot is None:  
                if  maxpsi/minpsi > 10.0 :
                    print("只有正的")
                    self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[0.1*limpsi], transparent=False ,vmax=0.1*limpsi)
    		    
                elif  maxpsi/minpsi < 0.1 :
                    print("只有负的")
                    self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[-0.1*limpsi], transparent=False ,vmin=-0.1*limpsi)
                else:
                    print("有正有负")
                    self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[-0.1*limpsi,0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
    		    
    
        #当场景被激活，或者参数发生改变，更新图像
        @on_trait_change(['posx','posy','posz','coef1','coef2'])
        def update_plot(self):
    	
            #print(self.coef1)
            #print(self.coef1)
    		#计算第二个波函数数值，把位置移到第二个波函数的原子核位置
            x2 = x-self.posx * a * n2**1.65
            y2 = y-self.posy * a * n2**1.65
            z2 = z-self.posz * a * n2**1.65
            r,Theta,Phi = cart2sph(x2,y2,z2)
    		
            scalars2 = calc_psi(r,Theta,Phi,n2,l2,m2)
    		
            scalars=scalars1*self.coef1 + scalars2*self.coef2
    
            maxpsi=np.abs(np.max(scalars))
            minpsi=np.abs(np.min(scalars))
            print("maxpsi",np.max(scalars))
            print("minpsi",np.min(scalars))
            limpsi=np.maximum(maxpsi,minpsi)
    
            mlab.clf()
            if  maxpsi/minpsi > 10.0 :
                print("只有正的")
                self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[0.1*limpsi], transparent=False ,vmax=0.1*limpsi)
    	
            elif  maxpsi/minpsi < 0.1 :
                print("只有负的")
                self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[-0.1*limpsi], transparent=False ,vmin=-0.1*limpsi)
            else:
                print("有正有负")
                self.plot = self.scene.mlab.contour3d(x,y,z,scalars, contours=[-0.1*limpsi,0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
    
    
    
        #建立视图布局
        view = View(
            Item("scene",editor=SceneEditor(scene_class=MayaviScene),
                 height=400,width=500,show_label=False),
            Group("_","posx","posy","posz","coef1","coef2"),
            resizable=True
        )
    
    
    
    print("请尽量保证两个波函数的主量子数接近")
    print("否则能量差别较大，无法有效重叠")
    print("例如")
    print("s-s轨道组合，选择1s-1s或1s-2s，不要选择1s-5s")
    print("s-p轨道组合，选择2s-2p或3s-3p，不要选择5s-2p")

    print("请输入第1个波函数的量子数")
    n1,l1,m1 = input_nlm()
    print("请输入第2个波函数的量子数")
    n2,l2,m2 = input_nlm()
    
    	
    
    #判断图像的范围
    #第一个波函数在坐标原点，第二在a*5*(n1^1.65+n2^1.65)的位置
    #范围选择要再大一点，left_side=-a*(5*n1**1.65),right_side=a*(5*n2**1.65)+a*5*(n1^1.65+n2^1.65)
    left_side=-a*(5*n1**1.65)
    right_side=a*(5*n2**1.65)+a*5*(n1**1.65+n2**1.65)
    
    #根据范围，确定格点
    x, y, z = np.ogrid[left_side:right_side:168j, left_side:right_side:168j, left_side:right_side:168j]
    r,Theta,Phi = cart2sph(x,y,z)
    x,y,z = sph2cart(r,Theta,Phi) #再次转换，保持xyz和psi有合适的维度
    scalars1 = calc_psi(r,Theta,Phi,n1,l1,m1)
    
    
    model = MyModel()
    model.configure_traits()

def Psi_iso():
    
    n,l,m = input_nlm()
    x, y, z = np.ogrid[-a*(5*n**1.65):a*(5*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j]
    r,Theta,Phi = cart2sph(x,y,z)
    x,y,z = sph2cart(r,Theta,Phi) #再次转换，保持xyz和psi有合适的维度

    
    scalars = calc_psi(r,Theta,Phi,n, l, m)
    #print('scalars')
    #print(scalars)
    
    maxpsi=np.abs(np.max(scalars))
    minpsi=np.abs(np.min(scalars))
    if n == 1:
        limpsi=maxpsi
    else:
        limpsi=np.minimum(maxpsi,minpsi)

#d轨道选择0.1有点小,可选择0.5,轨道更好看,更接近我们通常看到的轨道
    if n==1 and l== 0 :
        mlab.contour3d(x,y,z,scalars, contours=[0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
    else:
        mlab.contour3d(x,y,z,scalars, contours=[-0.1*limpsi,0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
#    mlab.volume_slice(scalars, plane_orientation='x_axes')#似乎可以用这个切面，待完成
#    mlab.colorbar(nb_colors=2)
#    mlab.scalarbar(nb_colors=2)

#PlotNodes
    mlab.contour3d(x, y, z, scalars, contours=[0], opacity=0.2, transparent=True, color=(0.5,0.5,0.0))
    
    #s轨道展示剖面，显示多层
    if l == 0:
        x, y, z = np.ogrid[0:a*(10*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j]
        r,Theta,Phi = cart2sph(x,y,z)
        x,y,z = sph2cart(r,Theta,Phi)
        x = x + a*(10*n**1.65)
        scalars = calc_psi(r,Theta,Phi,n, l, m)
        
        maxpsi=np.abs(np.max(scalars))
        minpsi=np.abs(np.min(scalars))
        if n == 1:
            limpsi=maxpsi
        else:
            limpsi=np.minimum(maxpsi,minpsi)
        
        if n == 1:
            mlab.contour3d(x-a*(10*n**1.65),y+a*(10*n**1.65),z,scalars, contours=[0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
        else:
            mlab.contour3d(x-a*(10*n**1.65),y+a*(10*n**1.65),z,scalars, contours=[-0.1*limpsi,0.1*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)

    mlab.title('isosurface of Psi')
    if l == 0:
        mlab.title('isosurface of Psi, semisphere shows the internal structure')
    mlab.show()


def Psi_iso_n():

    n,l,m = input_nlm()
    x, y, z = np.ogrid[-a*(5*n**1.65):a*(5*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j, -a*(5*n**1.65):a*(5*n**1.65):168j]
    r,Theta,Phi = cart2sph(x,y,z)
    x,y,z = sph2cart(r,Theta,Phi) #再次转换，保持xyz和psi有合适的维度
    
    
    scalars = calc_psi(r,Theta,Phi,n, l, m)
    #print('scalars')
    #print(scalars)
    
    maxpsi=np.abs(np.max(scalars))
    minpsi=np.abs(np.min(scalars))
    if n == 1:
        limpsi=maxpsi
    else:
        limpsi=np.minimum(maxpsi,minpsi)

    for nth_psi in range(1,6):
        if n==1 and l== 0 :
            mlab.contour3d(x,y+a*(10*nth_psi*n**1.65),z,scalars, contours=[0.1*nth_psi*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
        else:
            mlab.contour3d(x,y+a*(10*nth_psi*n**1.65),z,scalars, contours=[-0.1*nth_psi*limpsi,0.1*nth_psi*limpsi], transparent=False ,vmax=0.1*limpsi,vmin=-0.1*limpsi)
    #    mlab.colorbar(nb_colors=2)
    #    mlab.scalarbar(nb_colors=2)
    
    
    mlab.title('isosurface of Psi')
    mlab.show()


def Psi_prof():

    n,l,m = input_nlm()
    print("Select the projection plane")
    print("0: Z=0")
    print("1: X=0")
    print("2: Y=0")
    print("3: all")

    while True :
        try:
            planexyz = int(input())
            if (planexyz!=0) and (planexyz!=1) and (planexyz!=2) and (planexyz!=3):
                print('输入0,1,2,3分别代表的平面')
                continue
            else:
                break
        except:
            print('只能输入代表平面的整数')
    
    r = np.linspace(0.0, a*(5*n**1.65), 181)
    
    if planexyz==0 :
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)

        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,y,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,y,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red','blue'])
        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("yaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red','blue'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        
    if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
    if planexyz==1 :
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
        
    if planexyz==2 :
        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    
    
    if planexyz==3:
        ax1=plt.subplot(2,3,1,aspect='equal')
        ax2=plt.subplot(2,3,2,aspect='equal')
        ax3=plt.subplot(2,3,3,aspect='equal')
#    if planexyz==0 :
        r = np.linspace(0.0, a*(5*n**1.65), 181)
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        print(limpsi)
        #ax1.contour(x,y,f,levels=[0.0], colors=['black'])
        ax1.contour(x,y,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        
        ax1.set_xlabel("xaxis")#x轴上的名字
        ax1.set_ylabel("yaxis")#y轴上的名字
        ax1.set_xticks([])#关闭x刻度
        ax1.set_yticks([])#关闭y刻度
        ax1.set_title('profile of Psi_xy', fontsize = 14, fontweight = 'bold')
            
            
#        if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        ax2.set_xlabel("yaxis")#x轴上的名字
        ax2.set_ylabel("zaxis")#y轴上的名字
        ax2.set_xticks([])#关闭x刻度
        ax2.set_yticks([])#关闭y刻度
        ax2.set_title('profile of Psi_yz', fontsize = 14, fontweight = 'bold')
#        if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f,levels=[-0.1*limpsi,0.1*limpsi], colors=['red', 'blue'])

#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        ax3.set_xlabel("xaxis")#x轴上的名字
        ax3.set_ylabel("zaxis")#y轴上的名字
        ax3.set_xticks([])#关闭x刻度
        ax3.set_yticks([])#关闭y刻度
        ax3.set_title('profile of Psi_xz', fontsize = 14, fontweight = 'bold')
    plt.show()



def Psi_prof_n():
    
    n,l,m = input_nlm()

    print("Select the projection plane")
    print("0: Z=0")
    print("1: X=0")
    print("2: Y=0")
    print("3: all")

    while True :
        try:
            planexyz = int(input())
            if (planexyz!=0) and (planexyz!=1) and (planexyz!=2) and (planexyz!=3):
                print('输入0,1,2,3分别代表的平面')
                continue
            else:
                break
        except:
            print('只能输入代表平面的整数')
    
    r = np.linspace(0.0, a*(5*n**1.65), 181)
    
    if planexyz==0 :
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)

        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,y,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,y,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])
        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("yaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        
    if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
    if planexyz==1 :
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
        
    if planexyz==2 :
        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    
    
    if planexyz==3:
        ax1=plt.subplot(2,3,1,aspect='equal')
        ax2=plt.subplot(2,3,2,aspect='equal')
        ax3=plt.subplot(2,3,3,aspect='equal')
#    if planexyz==0 :
        r = np.linspace(0.0, a*(5*n**1.65), 181)
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)

        #ax1.contour(x,y,f,levels=[0.0], colors=['black'])
        ax1.contour(x,y,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        
        ax1.set_xlabel("xaxis")#x轴上的名字
        ax1.set_ylabel("yaxis")#y轴上的名字
        ax1.set_xticks([])#关闭x刻度
        ax1.set_yticks([])#关闭y刻度
        ax1.set_title('profile of Psi_xy', fontsize = 14, fontweight = 'bold')
            
            
#        if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        ax2.set_xlabel("yaxis")#x轴上的名字
        ax2.set_ylabel("zaxis")#y轴上的名字
        ax2.set_xticks([])#关闭x刻度
        ax2.set_yticks([])#关闭y刻度
        ax2.set_title('profile of Psi_yz', fontsize = 14, fontweight = 'bold')
#        if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f,levels=[-0.8*limpsi,-0.4*limpsi,-0.2*limpsi,-0.1*limpsi,-0.05*limpsi,-0.01*limpsi,0.01*limpsi,0.05*limpsi,0.1*limpsi,0.2*limpsi,0.4*limpsi,0.8*limpsi],colors = ['black'])

#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        ax3.set_xlabel("xaxis")#x轴上的名字
        ax3.set_ylabel("zaxis")#y轴上的名字
        ax3.set_xticks([])#关闭x刻度
        ax3.set_yticks([])#关闭y刻度
        ax3.set_title('profile of Psi_xz', fontsize = 14, fontweight = 'bold')
    plt.show()


def Psi_sq_prof():

        
    n,l,m = input_nlm()
    print("Select the projection plane")
    print("0: Z=0")
    print("1: X=0")
    print("2: Y=0")
    print("3: all")

    while True :
        try:
            planexyz = int(input())
            if (planexyz!=0) and (planexyz!=1) and (planexyz!=2) and (planexyz!=3):
                print('输入0,1,2,3分别代表的平面')
                continue
            else:
                break
        except:
            print('只能输入代表平面的整数')
    
    r = np.linspace(0.0, a*(5*n**1.65), 181)
    
    if planexyz==0 :
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)

        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,y,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,y,f*f,levels=[0.1*limpsi], colors=['blue'])
        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("yaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f*f,levels=[0.1*limpsi], colors=['blue'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        
    if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
    if planexyz==1 :
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(y,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(y,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        plt.xlabel("yaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
        
    if planexyz==2 :
        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #C = plt.contour(x,z,f,levels=[0.0], colors=['black'])
        C = plt.contour(x,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        plt.xlabel("xaxis")#x轴上的名字
        plt.ylabel("zaxis")#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        plt.axis('equal')
        plt.title('profile of Psi', fontsize = 14, fontweight = 'bold')
    
    
    if planexyz==3:
        ax1=plt.subplot(2,3,1,aspect='equal')
        ax2=plt.subplot(2,3,2,aspect='equal')
        ax3=plt.subplot(2,3,3,aspect='equal')
#    if planexyz==0 :
        r = np.linspace(0.0, a*(5*n**1.65), 181)
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        Theta = theta
        R,Phi = np.meshgrid(r,phi)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)

        #ax1.contour(x,y,f,levels=[0.0], colors=['black'])
        ax1.contour(x,y,f*f,levels=[0.01*limpsi], colors=['blue'])

        
        ax1.set_xlabel("xaxis")#x轴上的名字
        ax1.set_ylabel("yaxis")#y轴上的名字
        ax1.set_xticks([])#关闭x刻度
        ax1.set_yticks([])#关闭y刻度
        ax1.set_title('profile of Psi_xy', fontsize = 14, fontweight = 'bold')
            
            
#        if planexyz==1 :
        phi = pi/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax2.contour(y,z,f,levels=[0.0], colors=['black'])
        ax2.contour(y,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        ax2.set_xlabel("yaxis")#x轴上的名字
        ax2.set_ylabel("zaxis")#y轴上的名字
        ax2.set_xticks([])#关闭x刻度
        ax2.set_yticks([])#关闭y刻度
        ax2.set_title('profile of Psi_yz', fontsize = 14, fontweight = 'bold')
#        if planexyz==2 :
        phi = 0
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f*f,levels=[0.01*limpsi], colors=['blue'])

        phi = pi
        theta = np.linspace(0,pi,181)
        Phi = phi
        R,Theta = np.meshgrid(r,theta)
        f = calc_psi_prof(R,Theta,Phi,n,l,m)
        x,y,z = sph2cart(abs(R),Theta,Phi)
        maxpsi=np.abs(np.max(f))
        minpsi=np.abs(np.min(f))
        limpsi=np.maximum(maxpsi,minpsi)
        #ax3.contour(x,z,f,levels=[0.0], colors=['black'])
        ax3.contour(x,z,f*f,levels=[0.01*limpsi], colors=['blue'])

#        plt.clabel(C,inline = True,fontsize = 10) #显示等高线数值
        ax3.set_xlabel("xaxis")#x轴上的名字
        ax3.set_ylabel("zaxis")#y轴上的名字
        ax3.set_xticks([])#关闭x刻度
        ax3.set_yticks([])#关闭y刻度
        ax3.set_title('profile of Psi_xz', fontsize = 14, fontweight = 'bold')
    plt.show()
    
    
def sph_harm_sq_prof():


    l,m = input_lm()
    print("Select the projection plane")
    print("0: Z=0")
    print("1: X=0")
    print("2: Y=0")
    print("3: all")

    while True :
        try:
            planexyz = int(input())
            if (planexyz!=0) and (planexyz!=1) and (planexyz!=2) and (planexyz!=3):
                print('输入0,1,2,3分别代表的平面')
                continue
            else:
                break
        except:
            print('只能输入代表平面的整数')


    
    if planexyz==0:
        theta = pi/2
        phi = np.linspace(0,2*pi,181)

        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x1,y1,z1 = sph2cart(rholow*rholow,theta,phi)
        x2,y2,z2 = sph2cart(rhoup*rhoup,theta,phi)
        plt.axis('equal')
        plt.plot(x1, y1, color = 'red')
        plt.plot(x2, y2, color = 'blue')
    
    
    if planexyz==1 :#需要分phi = pi/2 or pi*3/2
        phi = pi/2
        theta = np.linspace(0,pi,181)
    
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(rholow*rholow,theta,phi)
        x4,y4,z4 = sph2cart(rhoup*rhoup,theta,phi)
        plt.axis('equal')
        plt.plot(y3, z3, color = 'red')
        plt.plot(y4, z4, color = 'blue')
        
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(rholow*rholow,theta,phi)
        x4,y4,z4 = sph2cart(rhoup*rhoup,theta,phi)
        plt.axis('equal')
        plt.plot(y3, z3, color = 'red')
        plt.plot(y4, z4, color = 'blue')
    
    
    if planexyz==2 :#需要分phi = 0 or pi
        phi = 0
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(rholow*rholow,theta,phi)
        x6,y6,z6 = sph2cart(rhoup*rhoup,theta,phi)
        plt.axis('equal')
        plt.plot(x5, z5, color = 'red')
        plt.plot(x6, z6, color = 'blue')
        
        phi = pi
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(rholow*rholow,theta,phi)
        x6,y6,z6 = sph2cart(rhoup*rhoup,theta,phi)
        plt.axis('equal')
        plt.plot(x5, z5, color = 'red')
        plt.plot(x6, z6, color = 'blue')



    if planexyz==0 :
        plt.xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Y_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    if planexyz==1 :
        plt.xlabel("Y_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    if planexyz==2 :
        plt.xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    
    if planexyz==3:
        ax1=plt.subplot(2,3,1,aspect='equal')
        ax2=plt.subplot(2,3,2,aspect='equal')
        ax3=plt.subplot(2,3,3,aspect='equal')
        
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x1,y1,z1 = sph2cart(rholow*rholow,theta,phi)
        x2,y2,z2 = sph2cart(rhoup*rhoup,theta,phi)
        ax1.axis('equal')
        ax1.plot(x1, y1, color = 'red')
        ax1.plot(x2, y2, color = 'blue')
    
    
    #if planexyz==1 or planexyz==3:#需要分phi = pi/2 or pi*3/2
        phi = pi/2
        theta = np.linspace(0,pi,181)
    
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(rholow*rholow,theta,phi)
        x4,y4,z4 = sph2cart(rhoup*rhoup,theta,phi)
        ax2.axis('equal')
        ax2.plot(y3, z3, color = 'red')
        ax2.plot(y4, z4, color = 'blue')
        
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(rholow*rholow,theta,phi)
        x4,y4,z4 = sph2cart(rhoup*rhoup,theta,phi)
        ax2.axis('equal')
        ax2.plot(y3, z3, color = 'red')
        ax2.plot(y4, z4, color = 'blue')
    
    
    #if planexyz==2 or planexyz==3:
        phi = 0
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(rholow*rholow,theta,phi)
        x6,y6,z6 = sph2cart(rhoup*rhoup,theta,phi)
        ax3.axis('equal')
        ax3.plot(x5, z5, color = 'red')
        ax3.plot(x6, z6, color = 'blue')
        
        phi = pi
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(rholow*rholow,theta,phi)
        x6,y6,z6 = sph2cart(rhoup*rhoup,theta,phi)
        ax3.axis('equal')
        ax3.plot(x5, z5, color = 'red')
        ax3.plot(x6, z6, color = 'blue')


    #if planexyz==0 or planexyz==3:
        ax1.set_xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax1.set_ylabel("Y_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax1.set_xticks([])#关闭x刻度
        ax1.set_yticks([])#关闭y刻度
        ax1.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    #if planexyz==1 or planexyz==3:
        ax2.set_xlabel("Y_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax2.set_ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax2.set_xticks([])#关闭x刻度
        ax2.set_yticks([])#关闭y刻度
        ax2.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    #if planexyz==2 or planexyz==3:
        ax3.set_xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax3.set_ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax3.set_xticks([])#关闭x刻度
        ax3.set_yticks([])#关闭y刻度
        ax3.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    
    plt.show()

def sph_harm_prof():

    l,m = input_lm()

    print("Select the projection plane")
    print("0: Z=0")
    print("1: X=0")
    print("2: Y=0")
    print("3: all")

    while True :
        try:
            planexyz = int(input())
            if (planexyz!=0) and (planexyz!=1) and (planexyz!=2) and (planexyz!=3):
                print('输入0,1,2,3分别代表的平面')
                continue
            else:
                break
        except:
            print('只能输入代表平面的整数')


    
    if planexyz==0:
        theta = pi/2
        phi = np.linspace(0,2*pi,181)

        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x1,y1,z1 = sph2cart(abs(rholow),theta,phi)
        x2,y2,z2 = sph2cart(abs(rhoup),theta,phi)
        plt.axis('equal')
        plt.plot(x1, y1, color = 'red')
        plt.plot(x2, y2, color = 'blue')
    
    
    if planexyz==1 :#需要分phi = pi/2 or pi*3/2
        phi = pi/2
        theta = np.linspace(0,pi,181)
    
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(abs(rholow),theta,phi)
        x4,y4,z4 = sph2cart(abs(rhoup),theta,phi)
        plt.axis('equal')
        plt.plot(y3, z3, color = 'red')
        plt.plot(y4, z4, color = 'blue')
        
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(abs(rholow),theta,phi)
        x4,y4,z4 = sph2cart(abs(rhoup),theta,phi)
        plt.axis('equal')
        plt.plot(y3, z3, color = 'red')
        plt.plot(y4, z4, color = 'blue')
    
    
    if planexyz==2 :#需要分phi = 0 or pi
        phi = 0
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(abs(rholow),theta,phi)
        x6,y6,z6 = sph2cart(abs(rhoup),theta,phi)
        plt.axis('equal')
        plt.plot(x5, z5, color = 'red')
        plt.plot(x6, z6, color = 'blue')
        
        phi = pi
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(abs(rholow),theta,phi)
        x6,y6,z6 = sph2cart(abs(rhoup),theta,phi)
        plt.axis('equal')
        plt.plot(x5, z5, color = 'red')
        plt.plot(x6, z6, color = 'blue')



    if planexyz==0 :
        plt.xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Y_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    if planexyz==1 :
        plt.xlabel("Y_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    if planexyz==2 :
        plt.xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        plt.ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        plt.xticks([])#关闭x刻度
        plt.yticks([])#关闭y刻度
        plt.title('profile of Y', fontsize = 14, fontweight = 'bold')
    
    if planexyz==3:
        ax1=plt.subplot(2,3,1,aspect='equal')
        ax2=plt.subplot(2,3,2,aspect='equal')
        ax3=plt.subplot(2,3,3,aspect='equal')
        
        theta = pi/2
        phi = np.linspace(0,2*pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x1,y1,z1 = sph2cart(abs(rholow),theta,phi)
        x2,y2,z2 = sph2cart(abs(rhoup),theta,phi)
        ax1.axis('equal')
        ax1.plot(x1, y1, color = 'red')
        ax1.plot(x2, y2, color = 'blue')
    
    
    #if planexyz==1 or planexyz==3:#需要分phi = pi/2 or pi*3/2
        phi = pi/2
        theta = np.linspace(0,pi,181)
    
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(abs(rholow),theta,phi)
        x4,y4,z4 = sph2cart(abs(rhoup),theta,phi)
        ax2.axis('equal')
        ax2.plot(y3, z3, color = 'red')
        ax2.plot(y4, z4, color = 'blue')
        
        phi = pi*3/2
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x3,y3,z3 = sph2cart(abs(rholow),theta,phi)
        x4,y4,z4 = sph2cart(abs(rhoup),theta,phi)
        ax2.axis('equal')
        ax2.plot(y3, z3, color = 'red')
        ax2.plot(y4, z4, color = 'blue')
    
    
    #if planexyz==2 or planexyz==3:
        phi = 0
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(abs(rholow),theta,phi)
        x6,y6,z6 = sph2cart(abs(rhoup),theta,phi)
        ax3.axis('equal')
        ax3.plot(x5, z5, color = 'red')
        ax3.plot(x6, z6, color = 'blue')
        
        phi = pi
        theta = np.linspace(0,pi,181)
        rho = calc_Y(theta,phi,l,m)
        rholow = np.ma.masked_where(rho < 0.0, rho)
        rhoup = np.ma.masked_where(rho >= 0.0, rho)
        x5,y5,z5 = sph2cart(abs(rholow),theta,phi)
        x6,y6,z6 = sph2cart(abs(rhoup),theta,phi)
        ax3.axis('equal')
        ax3.plot(x5, z5, color = 'red')
        ax3.plot(x6, z6, color = 'blue')


    #if planexyz==0 or planexyz==3:
        ax1.set_xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax1.set_ylabel("Y_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax1.set_xticks([])#关闭x刻度
        ax1.set_yticks([])#关闭y刻度
        ax1.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    #if planexyz==1 or planexyz==3:
        ax2.set_xlabel("Y_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax2.set_ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax2.set_xticks([])#关闭x刻度
        ax2.set_yticks([])#关闭y刻度
        ax2.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    #if planexyz==2 or planexyz==3:
        ax3.set_xlabel("X_axis",fontsize = 14, fontweight = 'bold')#x轴上的名字
        ax3.set_ylabel("Z_axis",fontsize = 14, fontweight = 'bold')#y轴上的名字
        ax3.set_xticks([])#关闭x刻度
        ax3.set_yticks([])#关闭y刻度
        ax3.set_title('profile of Y', fontsize = 14, fontweight = 'bold')
    
    plt.show()

def sph_harm_sq():
    ngrid = 131
    phi, theta = np.mgrid[0:2*pi:131j, 0:pi:131j]
    l,m = input_lm()
    s = calc_Y(theta,phi,l, m)
    x,y,z = sph2cart(s*s,theta,phi)
    
    
    # Represent spherical harmonics on the surface of the sphere
    
    if l == 0 :
        arr = np.ones((ngrid,ngrid), dtype=float)
        s = arr*s
    
    mlab.mesh(x , y , z, scalars=s/abs(s))
    mlab.title('Y^2,sph_harm^2')
    
    mlab.show()

def sph_harm():
    ngrid = 131
    l,m = input_lm()
    phi, theta = np.mgrid[0:2*pi:131j, 0:pi:131j]
    s = calc_Y(theta,phi,l, m)
    x,y,z = sph2cart(abs(s),theta,phi)
    
    
    # Represent spherical harmonics on the surface of the sphere
    
    if l == 0 :
        arr = np.ones((ngrid,ngrid), dtype=float)
        s = arr*s
    
    mlab.mesh(x , y , z, scalars=s/abs(s))
    mlab.title('Y,sph_harm')
    
    mlab.show()


def R_D():
    
    ngrid = 100
    n,l = input_nl()
    
    def plot_R():
        fig = plt.figure()
    
        #R = np.sqrt(X**2 + Y**2 + Z**2)
        R = np.linspace(0.0, a*(5*n**1.65), n*ngrid)
        rho = 2.0*R/n/a
        C = np.sqrt((2.0/n/a)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3))
        f = np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
        f = f*C
        R = R/a
        ax1=plt.subplot(2,3,1)
        equation = r'{R_{(r)}}'
        plt.title('$%s$'  %(equation),loc='center',  fontsize = 14, fontweight = 'bold')
        plt.title('a',loc='left',y=1,  fontsize = 24, fontweight = 'bold')
        plt.plot(R,f,color = 'red')


    #    只显示两个坐标轴的方法
    #    ax = plt.gca()
    #    ax.spines["right"].set_color("none")
    #    ax.spines["top"].set_color("none")
    #    ax.spines["bottom"].set_position(("data", 0))
    #    ax.spines["left"].set_position(("data", 0))
        plt.xlim(left=0)#方框左侧为零点
        if n-l!=1:
            plt.axhline(y=0.0, c="blue", ls="--", lw=1)#引入参考系
        else:
            plt.ylim(bottom=0)#方框底部为零点
        plt.axhline(y=0.0, c="blue", ls="--", lw=1)#引入参考系
        #plt.xticks([])
        plt.yticks([])
#        plt.show()
#        plt.close()
        
    def plot_R2():
#        fig = plt.figure()
    
        #R = np.sqrt(X**2 + Y**2 + Z**2)
        R = np.linspace(0.0, a*(5*n**1.65), n*ngrid)
        rho = 2.0*R/n/a
        C = np.sqrt((2.0/n/a)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3))
        f = np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
        f = f*C
        f = f**2
        f_max = np.max(f)
        f = n*f/f_max
        R = R/a
        ax2=plt.subplot(2,3,2)
        equation = r'R_{(r)}^2'
        plt.title('$%s$'  %(equation),loc='center',  fontsize = 14, fontweight = 'bold')
        plt.title('b',loc='left',y=1,  fontsize = 24, fontweight = 'bold')
        
        #节点的位置需要具体的x数值,等后续添加
        #min_indx=np.argmin(f)/ngrid*10*a#min value index2value
        #max_indx=0.5*np.max(f)
        #plt.annotate('节点',xy=(min_indx,0.0),xytext=(min_indx,max_indx),arrowprops=dict(arrowstyle="simple"))
        
        plt.plot(R,f,color = 'red')
        #plt.axhline(y=0.0, c="blue", ls="--", lw=1)#引入参考系
    
    #    只显示两个坐标轴的方法
    #    ax = plt.gca()
    #    ax.spines["right"].set_color("none")
    #    ax.spines["top"].set_color("none")
    #    ax.spines["bottom"].set_position(("data", 0))
    #    ax.spines["left"].set_position(("data", 0))
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        if l==0 and n!=1 :
            plt.ylim(0,0.2)
        #plt.xticks([])
        plt.yticks([])
#        plt.show()
#        plt.close()
        
    def plot_R3():
#        fig = plt.figure()
    
        #R = np.sqrt(X**2 + Y**2 + Z**2)
        R = np.linspace(0.0, a*(5*n**1.65), n*ngrid)
        rho = 2.0*R/n/a
        C = np.sqrt((2.0/n/a)**3*factorial(n-l-1)/(2*n*factorial(n+l)**3))
        f = np.exp(-rho/2)*rho**l*associatedLaguerre(l, n, rho)
        f = f*C
        f = R**2 * f**2
        R = R/a
        ax3=plt.subplot(2,3,3)
        equation = r'{r^2}R_{(r)}^2'
        plt.title('$%s$'  %(equation),loc='center',  fontsize = 14, fontweight = 'bold')
        plt.title('c',loc='left',y=1,  fontsize = 24, fontweight = 'bold')
        
        #节点的位置需要具体的x数值,等后续添加
        plt.plot(R,f,color = 'red')
        #plt.axhline(y=0.0, c="blue", ls="--", lw=1)#引入参考系
    
    #    只显示两个坐标轴的方法
    #    ax = plt.gca()
    #    ax.spines["right"].set_color("none")
    #    ax.spines["top"].set_color("none")
    #    ax.spines["bottom"].set_position(("data", 0))
    #    ax.spines["left"].set_position(("data", 0))
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        #plt.xticks([])
        plt.yticks([])

    plot_R()
    plot_R2()
    plot_R3()
#空白图片,图片中加入文本说明
    plt.subplot(2,3,4)
#loc只能这三个位置才会‘同时’出现三个title,否则取最后一个
#    plt.title("right bottom",y=0,loc='right', fontsize = 14, fontweight = 'bold')
#    plt.title("left top",y=0.8,loc='left', fontsize = 14, fontweight = 'bold')
#    plt.title("center",y=0.4,loc='center', fontsize = 14, fontweight = 'bold')
    if n==1 and l==0:
        equation = r'{R_{1,0}} = 2{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}\exp ( - \frac{{Zr}}{{{a_0}}})'
        line1='当n=1,l=0时,R函数的方程如下方公式所示\n图a,径向部分R反映了在任意角度方向上波函数随距离r变化的情况\n'
        line2='图b,径向部分的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,概率密度不为零,表明波函数中不出现径向节面\n'
        line3='图c,径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=1个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=0,代表s轨道,此时电子在核附近出现的概率密度不为零,这意味着处于s态的电子可以出现在原子核内,这也是Fermi接触作用的起源\n'
    elif n==2 and l==0:
        equation = r'{R_{2,0}} = \frac{1}{{2\sqrt {2} }}{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}\left( {2 - \frac{{Zr}}{{{a_0}}}} \right)\exp \left( { - \frac{{Zr}}{{2{a_0}}}} \right)'
        line1='当n=2,l=0时,R函数的方程如下方公式所示\n图a,径向部分R反映了在任意角度方向上波函数随距离r变化的情况；\n'
        line2='图b,径向部分的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,图中有1个概率密度为零的点,表明波函数中出现1个径向节面\n'
        line3='图c,径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=2个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=0,代表s轨道,此时电子在核附近出现的概率密度不为零,这意味着处于s态的电子可以出现在原子核内,这也是Fermi接触作用的起源\n'
    elif n==2 and l==1:
        equation = r'{R_{2,1}} = \frac{1}{{2\sqrt {6} }}{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}\left( {\frac{{Zr}}{{{a_0}}}} \right)\exp \left( { - \frac{{Zr}}{{2{a_0}}}} \right)'
        line1='当n=2,l=1时,R函数的方程如下方公式所示\n图a,径向部分R反映了在任意角度方向上波函数随距离r变化的情况\n'
        line2='图b,径向部分的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,图中无概率密度为零的点,表明波函数中不出现径向节面\n'
        line3='图c,径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=1个极大值峰,峰所在的位置就是电子出现概率密度极大值的位置\n量子数l=1,代表p轨道,此时电子在核附近出现的概率密度为零,这意味着处于p态的电子不可以出现在原子核内\n'
    elif n==3 and l==0:
        equation = r'{R_{3,0}} = \frac{2}{{81\sqrt {3} }}{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}\left( {27 - 18\frac{{Zr}}{{{a_0}}} + 2{{(\frac{{Zr}}{{{a_0}}})}^2}} \right)\exp ( - \frac{{Zr}}{{3{a_0}}})'
        line1='当n=3,l=0时,R函数的方程如下方公式所示\n图a:径向函数R反映了在任意给定的角度方向上波函数随r变化的情况\n'
        line2='图b:径向函数的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,图中有1个概率密度为零的点,表明波函数中出现1个径向节面\n'
        line3='图c:径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=3个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=0,代表s轨道,此时电子在核附近出现的概率密度不为零,这意味着处于s态的电子可以出现在原子核内,这也是Fermi接触作用的起源\n'
    elif n==3 and l==1:
        equation = r'{R_{3,1}} = \frac{4}{{81\sqrt {6} }}{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}\left( {6\frac{{Zr}}{{{a_0}}} - {{(\frac{{Zr}}{{{a_0}}})}^2}} \right)\exp ( - \frac{{Zr}}{{3{a_0}}})'
        line1='当n=3,l=1时,R函数的方程如下方公式所示\n图a:径向函数R反映了在任意给定的角度方向上波函数随r变化的情况\n'
        line2='图b:径向函数的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,图中有1个概率密度为零的点,表明波函数中出现1个径向节面\n'
        line3='图c:径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=2个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=1,代表p轨道,此时电子在核附近出现的概率密度为零,这意味着处于p态的电子不可以出现在原子核内\n'
    elif n==3 and l==2:
        equation = r'{R_{3,2}} = \frac{4}{{81\sqrt {30} }}{\left( {\frac{Z}{{{a_0}}}} \right)^{\frac{3}{2}}}{\left( {(\frac{{Zr}}{{{a_0}}})} \right)^2}\exp ( - \frac{{Zr}}{{3{a_0}}})'
        line1='当n=3,l=2时,R函数的方程如下方公式所示\n图a:径向函数R反映了在任意给定的角度方向上波函数随r变化的情况；\n'
        line2='图b:径向函数的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度,图中无概率密度为零的点,表明波函数中不出现径向节面\n'
        line3='图c:径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有n-l=1个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=2,代表d轨道,此时电子在核附近出现的概率密度为零,这意味着处于d态的电子不可以出现在原子核内\n'
    else:
        equation = r'{R_{n,l}}(r) = N{e^{ - \frac{\rho }{2}}}{\rho ^l}L_{n + l}^{2l + 1}(\rho ),\rho  = \frac{{2Zr}}{{n{a_0}}},N =  - {\left[ {{{(\frac{{2Z}}{{n{a_0}}})}^3} \cdot \frac{{(n - l - 1)!}}{{2n{{\left[ {(n + l)!} \right]}^3}}}} \right]^{\frac{1}{2}}},L_{n + l}^{2l + 1}(\rho ) = \frac{{{d^{2l + 1}}}}{{d{\rho ^{2l + 1}}}}\left[ {{e^\rho }\frac{{{d^{n + l}}}}{{d{\rho ^{n + l}}}}({e^{ - \rho }}{\rho ^{n + l}})} \right]'
        line1='当主量子数为n,角量子数为l时,R函数的方程如下方公式所示\n图a:径向函数R反映了在任意给定的角度方向上波函数随r变化的情况\n'
        line2='图b:径向函数的平方'r'$R_{(r)}^2$''反映了在距核r处,电子出现的径向概率密度\n'
        line3='图c:径向分布函数D,反映了在厚度dr的球壳内找到电子的概率密度,图中有（n-l）个极大值峰,峰所在的位置就是电子出现概率密度最大的位置\n量子数l=0,1,2,3,4,分别代表s,p,d,f,g轨道,ns态的电子在核附近的概率密度不为零,其他态的电子在核附近的概率密度为零\n\n'
#y=0.4,这个数值需要根据输出内容进行调节
    plt.title(line1+line2+line3+'$%s$'  %(equation),y=0.4,loc='left',  fontsize = 14, fontweight = 'bold',color = 'blue',verticalalignment='bottom')

    plt.plot()
    ax = plt.gca()
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    plt.xticks([])
    plt.yticks([])
    
    plt.show()
    plt.close()

print("SDUT orbital viewer")
print("""
          #######  ##     ## 
         ##     ## ##     ## 
         ##     ## ##     ## 
         ##     ## ##     ## 
         ##     ##  ##   ##  
         ##     ##   ## ##   
          #######     ###    
""")
print("Core Developer:牟景林; mjlink@126.com")
print("Contributors: 林文静,周佳宇")
print()

TIPS_realwf="""
本程序只能展示实波函数,
波函数量子数l和m在本程序中与图像的主要对应关系如下:
l=1,m=1,对应p_x轨道
l=1,m=-1,对应p_y轨道
l=2,m=-2,对应d_{x^2-y^2}轨道
l=2,m=-1,对应d_{yz}轨道
l=2,m=1,对应d_{xz}轨道
l=2,m=2,对应d_{xy}轨道
"""

TIPS_isosurf="""
波函数的等值面绘制中需要选取一定的数值,
数值过大过小可能出现部分等值面不能显示的情况，
具体情况可参考功能10，查看一个波函数的多个等值面。
"""


TIPS_test="""
该功能目前仅供测试，
可能有诸多bug。
"""

while True :
    print("1:R(r),R^2(r),D=r^2R^2(r),径向波函数,径向密度函数,径向分布函数")
    print("2:Y,Spherical harmonics,波函数角度部分/球谐函数Y")
    print("3:profile of Y,Y的剖面图")
    print("4:Y^2,Y的平方")
    print("5:profile of Y^2,Y平方的剖面图")
    print("6:contour of Psi,波函数的等值线图")
    print("7:contours of Psi,波函数的多个等值线")
    print("8:contour of Psi^2,概率密度的等值线图")
    print("9:isosurface of Psi,波函数的等值面")
    print("10:isosurfaces of Psi,波函数的多个等值面")
    print("11:probability density cloud,电子云")
    print("12:superposition of Psi,波函数的叠加")
    print("0:exit,退出")
    print("Please select the number to plot,请根据绘制的图形输入数字")
    try:
        plottype = int(input())
        if plottype==1 :
            R_D()
        if plottype==2 :
            print(TIPS_realwf)
            print()
            sph_harm()
        if plottype==3 :
            print(TIPS_realwf)
            print()
            sph_harm_prof()
            print()
        if plottype==4 :
            print(TIPS_realwf)
            print()
            sph_harm_sq()
            print()
        if plottype==5 :
            print(TIPS_realwf)
            print()
            sph_harm_sq_prof()
            print()
        if plottype==6 :
            print(TIPS_realwf)
            print()
            Psi_prof()
        if plottype==7 :
            print(TIPS_realwf)
            print()
            Psi_prof_n()
        if plottype==8 :
            print(TIPS_realwf)
            print()
            print(TIPS_isosurf)
            print()
            Psi_sq_prof()
        if plottype==9 :
            print(TIPS_realwf)
            print()
            print(TIPS_isosurf)
            print()
            Psi_iso()
        if plottype==10 :
            print(TIPS_realwf)
            print()
            Psi_iso_n()
        if plottype==11 :
            print(TIPS_realwf)
            print()
            Psi_points()
        if plottype==12 :
            print(TIPS_realwf)
            print()
            print(TIPS_isosurf)
            Psi_superpos()
        if plottype==0 :
            break
        else:
            continue
    except ValueError:
        print("只能输入对应的数字")


