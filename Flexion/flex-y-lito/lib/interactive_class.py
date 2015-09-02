import numpy as np
from flex_y_lito_functions import *
from double_window_selection import DoubleWindowSelection
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText


class InteractivePlot(object):
    def __init__(self, x, y, topo, w, t, padding, (rho_t,rho_c,rho_m,E,nu),
                 (Te_min,Te_max,n_Te), (topo_ax,w_ax,rms_ax,win_ax1,win_ax2),
                  outfilename):
        self.x = x
        self.y = y
        self.topo = topo
        self.w = w
        self.t = t
        self.padding = padding
        self.rho_t = rho_t
        self.rho_c = rho_c
        self.rho_m = rho_m
        self.E, self.nu = E, nu
        self.Te_min, self.Te_max, self.n_Te = Te_min, Te_max, n_Te
        self.topo_ax = topo_ax
        self.w_ax = w_ax
        self.rms_ax = rms_ax
        self.win_ax1, self.win_ax2 = win_ax1, win_ax2
        self.vmin, self.vmax = 0, np.max(self.w)
        self.outfilename = outfilename
        
        print "\nFLEXURE INSTRUCTIONS:"
        print "Press Enter or Right Click to calculate the elastic thickness that fits the window data."
        print "Press t to change the Te axe (Te_min, Te_max, ammount of Te points)."
        print "Press v to adjust the min/max values of the windows' colormaps"
        print "Press z to save the fitting data.\n"
        
        
        
        self.Te_fit, self.rms = None, None
        
        self.win_ax1.set_xticklabels([])
        self.win_ax1.set_yticklabels([])
        self.win_ax2.set_xticklabels([])
        self.win_ax2.set_yticklabels([])
        self.win_ax1.set_xlabel("Ventana de w Invertida")
        self.win_ax2.set_xlabel("Ventana de w por Flexion")
        
        self.x_done = []
        self.y_done = []
        self.l_done_topo, = self.topo_ax.plot(None, None, 'o', color='k', ms=2)
        self.l_done_w, = self.w_ax.plot(None, None, 'o', color='k', ms=2)
        
        
        self.l, = self.rms_ax.plot(None,None,'.')
        self.rms_ax.set_xlabel("Te")
        self.rms_ax.set_ylabel("RMS")
        self.ny, self.nx = np.shape(x)
        self.dx = (np.max(self.x)-np.min(self.x))/(self.nx-1)
        self.window = DoubleWindowSelection(self.x,self.y,self.topo_ax,self.w_ax)
        self.topo_ax.figure.canvas.mpl_connect('key_press_event', self.on_key)
        self.topo_ax.figure.canvas.mpl_connect('button_press_event', self.on_press)
        
        
    def on_key(self, event):
        if event.key == "enter":
            self.Te_calculation_on_plot()
        if event.key == "t":
            h = raw_input("Te min [" + str(self.Te_min) + "]: ")
            if h != "":
                self.Te_min = float(h)
            h = raw_input("Te max ["+ str(self.Te_max) + "]: ")
            if h != "":
                self.Te_max = float(h)
            h = raw_input("Ammount of Te points [" + str(self.n_Te) + "]: ")
            if h != "":
                self.n_Te = int(h)
            self.Te_calculation_on_plot()
        if event.key == "z":
            self.Te_calculation_on_plot()
            self.save_data()
        if event.key == "v":
            h = raw_input("Min value of w to plot in the windows: [" + str(self.vmin) + "]: ")
            if h != "":
                self.vmin = float(h)
            h = raw_input("Max value of w to plot in the windows: [" + str(self.vmax) + "]: ")
            if h != "":
                self.vmax = float(h)
            self.Te_calculation_on_plot()
    
    def on_press(self, event):
        if event.button == 3:
            self.Te_calculation_on_plot()
    
    def save_data(self):
        ## Write file
        line = str(self.window.x_center) + "\t" + str(self.window.y_center) + "\t" + \
                str(2*self.window.half_width) + "\t" + str(self.Te_fit) + "\t" + \
                str(self.rms) + "\n"
        print "\nx_center\ty_center\twin width\tTe\tw rms"
        print line
        outfile = open(self.outfilename,'a')
        outfile.write(line)
        outfile.close()
        ## Refresh dots
        self.x_done.append(self.window.x_center)
        self.y_done.append(self.window.y_center)
        self.refresh_done_pts()
        
    
    def Te_calculation_on_plot(self):
        ## Get the cutted inverted deflexion
        area = (self.window.x1, self.window.x2, self.window.y1, self.window.y2)
        x_win, y_win, [w_inverted_cut] = inside_points(self.x, self.y, [self.w], area)

        ## Get posibles roots
        Te_list, w_list = w_possible_list(self.x, self.y, self.topo,
                                          self.padding, self.t,
                                          self.rho_t, self.rho_c, self.rho_m,
                                          self.E, self.nu,
                                          self.Te_min, self.Te_max, self.n_Te
                                          )
        ## Cut the deflexions
        x_win, y_win, w_list_cut = inside_points(self.x, self.y, w_list, area)
        
        ## Compare with inverted deflexion
        rms = [rms_calculation(w_inverted_cut, w_i) for w_i in w_list_cut]
        Te_fit = Te_list[rms.index(min(rms))]
        w_fit_cut = w_list_cut[rms.index(min(rms))]
        
        self.Te_fit = Te_fit
        self.rms = np.min(rms)
        
        ## Update plot
        self.win_ax1.cla()
        self.win_ax1.contourf(x_win, y_win, w_inverted_cut, 150, vmin=self.vmin, vmax=self.vmax)
        self.win_ax1.set_xticklabels([])
        self.win_ax1.set_yticklabels([])
        self.win_ax2.cla()
        self.win_ax2.contourf(x_win, y_win, w_fit_cut, 150, vmin=self.vmin, vmax=self.vmax)
        self.win_ax2.set_xticklabels([])
        self.win_ax2.set_yticklabels([])
        self.win_ax1.set_xlabel("Ventana de w Invertida")
        self.win_ax2.set_xlabel("Ventana de w por Flexion")
        

        
        
        self.l.set_xdata(Te_list)
        self.l.set_ydata(rms)
        self.rms_ax.set_title("Te="+str(Te_fit))
        self.rms_ax.set_xlim(self.Te_min, self.Te_max)
        self.rms_ax.set_ylim(min(rms)*0.8, max(rms)*1.2)
        self.rms_ax.figure.canvas.draw()


    def refresh_done_pts(self):
        self.l_done_topo.set_xdata(self.x_done)
        self.l_done_topo.set_ydata(self.y_done)
        self.l_done_w.set_xdata(self.x_done)
        self.l_done_w.set_ydata(self.y_done)
        self.topo_ax.figure.canvas.draw()
        self.w_ax.figure.canvas.draw()
        
        
