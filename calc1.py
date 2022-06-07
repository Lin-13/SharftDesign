from turtle import shape
from cv2 import idct
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt
from calc import theta,deg,a,b,e,y1,y2
import calc
x0=calc.y
phi=[deg[np.min(calc.x2)],deg[np.max(calc.x2)],\
    deg[np.min(calc.x1)],deg[np.max(calc.x1)]]
x0_min=np.min(x0)
print(phi)
L_left=30#滑块到左夹板距离
x1=x0+L_left#左边夹板轨迹
L=20#弹簧长度
x2=np.zeros(x0.shape)
deltaL_base=-1
deltaL=3*deltaL_base#弹簧形变，2倍余量
a=0.9
delta=3
#右夹板轨迹
for i in range(len(x0)):
    x2[i]=x0[i]+L_left+2*L+50
    if x0[i]<y2:#离开deltaL逐渐变为0,delta变大
        x2[i]=a*y2+(1-a)*x0[i]+L_left+2*L+50+deltaL*(np.max(calc.x2)-i)/(np.max(calc.x2)-np.min(calc.x2))\
            +delta*(np.min(calc.x2)-i)/(np.min(calc.x2)-np.max(calc.x2))
    if x0[i]>y1:#夹紧delta->0,deltaL偏离
        x2[i]=a*y1+(1-a)*x0[i]+L_left+2*L+50+deltaL*(i-np.min(calc.x1))/(np.max(calc.x1)-np.min(calc.x1))\
            +delta*(np.max(calc.x1)-i)/(np.max(calc.x1)-np.min(calc.x1))
    if i>np.max(calc.x2) and i<np.min(calc.x1):#处于回程,delta变小
        x2[i]=x2[i]+delta
    if i<np.min(calc.x2) or i>np.max(calc.x1):#处于推程
        x2[i]=x2[i]+deltaL
            
#光滑轨迹
pointNum=10#采样数
f1=CubicSpline(theta,x2)#原曲线
theta2=np.linspace(0,np.pi*2,pointNum+1,endpoint=True)#采样点
rho2=f1(theta2)
rho2[-1]=rho2[0]
f2=CubicSpline(theta2,rho2,bc_type='periodic')#根据采样点插值函数
datapoint=20#插值数
theta_final=np.linspace(0,np.pi*2,datapoint+1,endpoint=True)
x2_final=f2(theta_final)#插值点
#确定凸轮形状

print("rho_final_min:",np.min(x2_final))
###
L_right=175#滑块中心到右夹板距离
###
baseR=np.min(x2_final)-L_right
print("baseR:",baseR)
rho_final=x2_final-L_right#导杆轨迹
x=rho_final*np.cos(theta_final)
y=rho_final*np.sin(theta_final)
R=1.5#凸轮轴承半径
x_new,y_new=x,y
for i in range(len(x)-1):
    k=(y[i+1]-y[i])/(x[i+1]-x[i])
    k_phi=np.arctan(k)
    if x[i+1]<x[i]:
        k_phi=np.pi+k_phi
    x_new[i]=x[i]+R*np.cos(k_phi+np.pi/2)
    y_new[i]=y[i]+R*np.sin(k_phi+np.pi/2)
x=x_new
y=y_new
x[-1]=x[0]
y[-1]=y[0]
#输出
file=open('data.txt','w')
#file.write('X Y Z\n')
for i in range(len(x)):
    file.write("{}mm {}mm {}mm\n".format(x[i],y[i],0))
#file.write("{}mm {}mm {}mm\n".format(x[0],y[0],0))
file.close()
###

###

x2_final_all=f2(theta)
plt.figure(2)#左右轨迹
plt.plot(deg,x1,'r')
#plt.plot(deg,x2,'r')
plt.plot(deg,x2_final_all,'b')
plt.legend(['left','right'])
plt.savefig("运动循环图.png")

x_space=x2-x1
x_space_final=f2(theta)-x1
plt.figure(3)#光滑
plt.title("space")
plt.plot(deg,x_space,'r')
plt.plot(deg,x_space_final,'b')

baseline=(50+2*L)*np.ones(len(theta))
baseline1=(50+2*L+deltaL_base)*np.ones(len(theta))
plt.plot(deg,baseline,'g')
plt.plot(deg,baseline1,'g')
plt.legend(['space','space_final','baseline','baseline_final'])
#plt.legend(["原函数","光滑后实际两板距离","接触基准","夹紧基准"])
plt.figure(4)#间距
plt.plot(deg,x_space_final,'b')

#轨迹参数
x2_min=np.min(x2_final_all)
x2_max=np.max(x2_final_all)
x2_min_where,x2_max_where=deg[np.argmin(x2_final_all)],deg[np.argmax(x2_final_all)]
x1_min,x1_max=np.min(x1),np.max(x1)
x1_min_where=deg[np.argmin(x1)]
x1_max_where=deg[np.argmax(x1)]
x_space_min=np.min(x_space_final)
x_space_max=np.max(x_space_final)
x_space_min_where,x_space_max_where=deg[np.argmin(x_space_final)],deg[np.argmax(x_space_final)]
print("x2_final_min:{},x2_final_max:{}".\
    format(x2_min,x2_max))
print("x2_min_where:{},x2_max_where:{}".\
    format(x2_min_where,x2_max_where))
print("x1_min:{},x1_max:{}".format(x1_min,x1_max))
print("x1_min_where:{},x1_max_where:{}".\
    format(x1_min_where,x1_max_where))
print("x_space_min:{},x_space_max:{}".format(x_space_min,x_space_max))
print("x_space_min_where:{},x_space_max_where:{}"\
    .format(x_space_min_where,x_space_max_where))
plt.figure(2)#画出最大最小值
plt.plot([x1_min_where,x1_min_where],[x1_min-1,x1_min],':')
plt.plot([deg[0],x1_min_where],[x1_min,x1_min],':')

plt.plot([x1_max_where,x1_max_where],[x1_min-1,x1_max],':')
plt.plot([deg[0],x1_max_where],[x1_max,x1_max],':')

plt.plot([x2_min_where,x2_min_where],[x1_min-1,x2_min],':')
plt.plot([deg[0],x2_min_where],[x2_min,x2_min],':')

plt.plot([x2_max_where,x2_max_where],[x1_min-1,x2_max],':')
plt.plot([deg[0],x2_max_where],[x2_max,x2_max],':')
plt.grid()
plt.savefig("左右轨迹.png")

plt.figure(4)#画出最大最小值-space
plt.plot([x_space_min_where,x_space_min_where],[x_space_min-1,x_space_min],':')
plt.plot([deg[0],x_space_min_where],[x_space_min,x_space_min],':')
plt.plot([x_space_max_where,x_space_max_where],[x_space_min-1,x_space_max],':')
plt.plot([deg[0],x_space_max_where],[x_space_max,x_space_max],':')
plt.grid()
plt.savefig("间距.png")
plt.figure(5)#各个阶段
plt.plot(deg,x_space_final,'b')
plt.grid()
plt.plot([deg[np.min(calc.x2)],deg[np.min(calc.x2)]],\
    [np.min(x_space_final)-1,x_space_final[np.min(calc.x2)]],'o')
plt.plot([deg[np.max(calc.x2)],deg[np.max(calc.x2)]],\
    [np.min(x_space_final)-1,x_space_final[np.max(calc.x2)]],'o')
plt.plot([deg[np.min(calc.x1)],deg[np.min(calc.x1)]],\
    [np.min(x_space_final)-1,x_space_final[np.min(calc.x1)]],'o')
plt.plot([deg[np.max(calc.x1)],deg[np.max(calc.x1)]],\
    [np.min(x_space_final)-1,x_space_final[np.max(calc.x1)]],'o')

plt.plot([deg[0],deg[np.min(calc.x2)]],[np.min(x_space_final)-1]*2,':or')
plt.plot([deg[np.min(calc.x2)],deg[np.max(calc.x2)]],\
    [np.min(x_space_final)-1]*2,':og')
plt.plot([deg[np.max(calc.x2)],deg[np.min(calc.x1)]],\
    [np.min(x_space_final)-1]*2,':ob')
plt.plot([deg[np.min(calc.x1)],deg[np.max(calc.x1)]],\
    [np.min(x_space_final)-1]*2,':oy')
plt.plot([deg[-1],deg[np.max(calc.x1)]],[np.min(x_space_final)-1]*2,':or')
#plt.legend(['1','2','3','4'])
plt.savefig("各个阶段.png")
plt.show()