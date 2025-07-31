import numpy as np
import pandas as pd
import time
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
mpl.rcParams.update({"font.family":"Times New Roman", "mathtext.fontset": "dejavuserif", "font.size": 16})


def plotter(t_shift_sec, t_tempMemory_s):
    fig, axes = plt.subplots(2, 2, figsize=(10, 5), sharex=True)
    [ax0, ax1], [ax2, ax3] = axes

    data = pd.DataFrame()

    def animate(i):
        global data

        for _ in range(10):
            # try to read the values for plotting 10 times
            # to prevent errors with opening the file
            try:
                df = pd.read_csv('buffer.txt',
                                 sep = r',',
                                 engine = 'python')
                break
            except:
                continue

        if i == 0:
            data = df.copy()
        else:
            data = data.append(df, ignore_index=True)

        if data.shape[0] > int(t_shift_sec / t_tempMemory_s):
            data.drop(labels=0, inplace=True)
            data.reset_index(drop=True, inplace=True)

        ax0.cla(), ax1.cla(), ax2.cla(), ax3.cla()
        
        # plot pressure
        line0, = ax0.plot(data['time/hh:mm:ss'],
                          data['PT101/MPa'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'k')
        line1, = ax2.plot(data['time/hh:mm:ss'],
                          data['DPT102/MPa'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'k')
        # plot Joule-Thomson temperature
        line2, = ax1.plot(data['time/hh:mm:ss'],
                          data['TT101/K'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'b',
                          label = 'TT101')
        line3, = ax1.plot(data['time/hh:mm:ss'],
                          data['TT102/K'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'r',
                          label = 'TT102')
        line4, = ax1.plot(data['time/hh:mm:ss'],
                          data['TT008/K'],
                          linestyle = '-.',
                          linewidth = 0.8,
                          color = 'g',
                          label = 'cold head')
        # plot HX temperature
        line5, = ax3.plot(data['time/hh:mm:ss'],
                          data['TT006/K'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'b',
                          label = 'HP in')
        line6, = ax3.plot(data['time/hh:mm:ss'],
                          data['TT007/K'],
                          linestyle = '-',
                          linewidth = 0.8,
                          color = 'r',
                          label = 'HP out')
        line7, = ax3.plot(data['time/hh:mm:ss'],
                          data['TT009/K'],
                          linestyle = '-.',
                          linewidth = 0.8,
                          color = 'b',
                          label = 'LP in')
        line8, = ax3.plot(data['time/hh:mm:ss'],
                          data['TT010/K'],
                          linestyle = '-.',
                          linewidth = 0.8,
                          color = 'r',
                          label = 'LP out')

        ax0.set_ylabel(r'$p~/~{\rm MPa}$')
        ax2.set_xlabel(r'$t~/~{\rm hh:mm:ss}$')
        ax2.set_ylabel(r'$\Delta p~/~{\rm MPa}$')
        
        ax1.set_ylabel(r'$T~/~{\rm K}$')
        ax3.set_xlabel(r'$t~/~{\rm hh:mm:ss}$')
        ax3.set_ylabel(r'$T~/~{\rm K}$')

        t = list(data['time/hh:mm:ss'])
        if len(t) > 2:
            ax0.set_xticks([t[0], t[int(len(t)/2)], t[-1]])
            ax1.set_xticks([t[0], t[int(len(t)/2)], t[-1]])
            ax2.set_xticks([t[0], t[int(len(t)/2)], t[-1]])
            ax3.set_xticks([t[0], t[int(len(t)/2)], t[-1]])

        ax1.legend(frameon=False)
        ax3.legend(frameon=False)
        
        plt.tight_layout(pad=0.2)
        return line0, line1, line2, line3, line4, line5, line6, line7, line8,

    anim = FuncAnimation(fig, animate, interval=1000*t_tempMemory_s)
    plt.show()


def main():
    t_shift_sec=180
    t_tempMemory_s=0.5
    plotter(t_shift_sec=t_shift_sec, t_tempMemory_s=t_tempMemory_s)


if __name__ == '__main__':
    main()

