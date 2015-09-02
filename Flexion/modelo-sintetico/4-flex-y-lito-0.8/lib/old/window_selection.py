import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def nearest_point(x,y,xp,yp):
    x_nearest, y_nearest = x[0][0], y[0][0]
    distance_x = abs(xp - x_nearest)
    distance_y = abs(yp - y_nearest)
    ny, nx = np.shape(x)
    for i in range(nx):
        if abs(x[0][i] - xp) < distance_x:
            x_nearest = x[0][i]
            distance_x = abs(xp - x_nearest)
    for j in range(ny):
        if abs(y[j][0] - yp) < distance_y:
            y_nearest = y[j][0]
            distance_y = abs(yp - y_nearest)
    return x_nearest, y_nearest
                

class WindowSelection(object):
    def __init__(self,x,y,ax):
        self.x = x
        self.y = y
        #~ self.z = z
        self.dx = abs(self.x[0][0] - self.x[0][-1])/(np.shape(self.x)[1]-1)
        self.dy = abs(self.y[0][0] - self.y[-1][0])/(np.shape(self.y)[0]-1)
        assert self.dx == self.dy, "dx != dy"
        self.rect = Rectangle((0,0), 1, 1,fc='None')
        self.x_center, self.y_center = None, None
        self.half_width = (min(np.shape(self.x))/16)*self.dx
        self.x1, self.x2 = None, None
        self.y1, self.y2 = None, None
        self.ax = ax
        self.l, = self.ax.plot([self.x_center],[self.y_center],'o')
        self.ax.add_patch(self.rect)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax.figure.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.ax.figure.canvas.mpl_connect('key_press_event', self.on_key)

    def on_press(self,event):
        if event.inaxes!=self.ax:
            return
        if event.button == 1:
            self.x_center, self.y_center = event.xdata, event.ydata
            self.x_center, self.y_center = nearest_point(self.x, self.y,
                                                  self.x_center, self.y_center)
            self.rectangle_construction()    
        
    def on_scroll(self,event):
        self.half_width += event.step * self.dx
        self.rectangle_construction()
        
    def on_key(self,event):
        event_list = ["right","left","up","down"]
        if event.key in event_list:
            if event.key == "right":
                self.x_center += self.dx
            elif event.key == "left":
                self.x_center -= self.dx
            elif event.key == "up":
                self.y_center += self.dx
            elif event.key == "down":
                self.y_center -= self.dx
            self.rectangle_construction()
        if event.key == "i":
            print "(x,y)=",(self.x_center,self.y_center),"Width:",self.half_width*2
        if event.key == "+" or event.key == "-":
            if event.key == "+":
                self.half_width += self.dx
            elif event.key == "-":
                self.half_width -= self.dx
            self.rectangle_construction()
        
    def rectangle_construction(self):
        self.x1 = self.x_center - self.half_width
        self.x2 = self.x_center + self.half_width
        self.y1 = self.y_center - self.half_width
        self.y2 = self.y_center + self.half_width
        self.rect.set_width(self.x2 - self.x1)
        self.rect.set_height(self.y2 - self.y1)
        self.rect.set_xy((self.x1, self.y1))
        self.l.set_xdata([self.x_center])
        self.l.set_ydata([self.y_center])
        self.ax.figure.canvas.draw()
