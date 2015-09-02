import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle              



class DoubleWindowSelection(object):
    def __init__(self,x,y,ax1,ax2):
        self.x = x
        self.y = y
        #~ self.z = z
        self.dx = abs(self.x[0][0] - self.x[0][-1])/(np.shape(self.x)[1]-1)
        self.dy = abs(self.y[0][0] - self.y[-1][0])/(np.shape(self.y)[0]-1)
        assert self.dx == self.dy, "dx != dy"
        self.rect1 = Rectangle((0,0), 1, 1,fc='None')
        self.rect2 = Rectangle((0,0), 1, 1,fc='None')
        self.x_center, self.y_center = None, None
        self.half_width = (min(np.shape(self.x))/8)*self.dx
        self.x1, self.x2 = None, None
        self.y1, self.y2 = None, None
        self.ax1 = ax1
        self.ax2 = ax2
        
        self.l1, = self.ax1.plot([self.x_center],[self.y_center],'o')
        self.l2, = self.ax2.plot([self.x_center],[self.y_center],'o')
        self.ax1.add_patch(self.rect1)
        self.ax2.add_patch(self.rect2)
        self.ax1.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.ax1.figure.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.ax1.figure.canvas.mpl_connect('key_press_event', self.on_key)
        
        print "\nWINDOWS INSTRUCTIONS:"
        print "Click to select the window center"
        print "Move the center with arrows or click again"
        print "Resize the window with the mouse scroll or with '+' and '-'"
        print "Press 'i' to show information about the window"
        

    def on_press(self,event):
        if event.inaxes == self.ax1 or event.inaxes == self.ax2:
            if event.button == 1:
                self.x_center, self.y_center = event.xdata, event.ydata
                #~ self.x_center, self.y_center = nearest_point(self.x, self.y,
                                                      #~ self.x_center, self.y_center)
                self.x_center -= self.x_center%self.dx
                self.y_center -= self.y_center%self.dx
                self.rectangle_construction()    
        else:
            return
            
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
        self.rect1.set_width(self.x2 - self.x1)
        self.rect1.set_height(self.y2 - self.y1)
        self.rect1.set_xy((self.x1, self.y1))
        self.l1.set_xdata([self.x_center])
        self.l1.set_ydata([self.y_center])
        
        self.rect2.set_width(self.x2 - self.x1)
        self.rect2.set_height(self.y2 - self.y1)
        self.rect2.set_xy((self.x1, self.y1))
        self.l2.set_xdata([self.x_center])
        self.l2.set_ydata([self.y_center])
        self.ax1.figure.canvas.draw()
 
