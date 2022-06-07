import numpy as np
import matplotlib.pyplot as plt

theta=np.linspace(np.pi/6,np.pi*2+np.pi/6,3600)
deg=theta/np.pi*180
a=30
b=100
e=30
#阈值
y1=122
y2=75
y=a*np.cos(theta)+np.sqrt(b**2-(e+a*np.sin(theta))**2)
y_max=np.max(y)
x_max=deg[np.where(y==y_max)]
print("y_max={},x={}".format(y_max,x_max))
y_min=np.min(y)
x_min=deg[np.where(y==y_min)]
print("y_min={},x={}".format(y_min,x_min))

x1=np.where(y>y1)
print("x1={},{}".format(deg[np.min(x1)],deg[np.max(x1)]))
x2=np.where(y<y2)
print("x2={},{}".format(deg[np.min(x2)],deg[np.max(x2)]))
plt.figure(1)
plt.plot(deg,y)
up=np.ones(len(deg))*y1
plt.plot(deg,up)
low=np.ones(len(deg))*y2
plt.plot(deg,low)
plt.xlabel('theta(°)')
plt.ylabel('y(mm)')
#plt.title("y=a*cos(theta)+sqrt(b**2-(e+a*sin(theta))**2),a={},b={},e={}".format(a,b,e))

plt.savefig("calc.png")
#plt.show()