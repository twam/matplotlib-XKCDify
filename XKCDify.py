# XKCD plot generator
# 
# Licensed under CC-BY-SA 3.0 (http://creativecommons.org/licenses/by-sa/3.0/)
#
# based on work by 
#   Jake Vanderplas (http://jakevdp.github.com/blog/2012/10/07/xkcd-style-plots-in-matplotlib)
#   Damon McDougall (http://www.mail-archive.com/matplotlib-users@lists.sourceforge.net/msg25499.html)

import numpy
import pylab
import matplotlib.font_manager
import scipy.interpolate
import scipy.signal

# We need a special font for the code below.  It can be downloaded this way:
import os
import urllib2
if not os.path.exists('Humor-Sans.ttf'):
    fhandle = urllib2.urlopen('http://antiyawn.com/uploads/Humor-Sans.ttf')
    open('Humor-Sans.ttf', 'wb').write(fhandle.read())

# Mimic a hand-drawn line from (x, y) data
def xkcd_line(ax, x, y, transform, mag=1.0, f1=30, f2=0.05, f3=15):
    x = numpy.asarray(x)
    y = numpy.asarray(y)

    x_scaled = numpy.copy(x)
    y_scaled = numpy.copy(y)

    # define transformation
    trans = ax.transScale+ax.transLimits

    # transform to axis coordinate systems
    if (transform != ax.transAxes) and (transform != 'none'):
        for i in range(len(x)):
            [[x_scaled[i], y_scaled[i]]] = trans.transform(numpy.array([(x[i], y[i])]))

    # compute the total distance along the path
    dx = x_scaled[1:] - x_scaled[:-1]
    dy = y_scaled[1:] - y_scaled[:-1]
    dist_tot = numpy.sum(numpy.sqrt(dx * dx + dy * dy))

    # number of interpolated points is proportional to the distance
    Nu = int(200 * dist_tot)
    u = numpy.arange(-1, Nu + 1) * 1. / (Nu - 1)

    # interpolate curve at sampled points
    k = min(3, len(x) - 1)
    res = scipy.interpolate.splprep([x_scaled, y_scaled], s=0, k=k)
    x_int, y_int = scipy.interpolate.splev(u, res[0]) 

    # we'll perturb perpendicular to the drawn line
    dx = x_int[2:] - x_int[:-2]
    dy = y_int[2:] - y_int[:-2]
    dist = numpy.sqrt(dx * dx + dy * dy)

    # create a filtered perturbation
    coeffs = mag * numpy.random.normal(0, 0.01, len(x_int) - 2)
    b = scipy.signal.firwin(f1, f2 * dist_tot, window=('kaiser', f3))
    response = scipy.signal.lfilter(b, 1, coeffs)

    x_int[1:-1] += response * dy / dist
    y_int[1:-1] += response * dx / dist

    # retransform to original coordinate system
    if (transform != ax.transAxes) and (transform != 'none'):
        for i in range(len(x_int)):
            [[x_int[i], y_int[i]]] = trans.inverted().transform(numpy.array([(x_int[i], y_int[i])]))
    
    return x_int[1:-1], y_int[1:-1]

# XKCDify plot
#
# Parameters
# ----------
# ax : Axes instance
#     the axes to be modified.
# mag : float
#     The magnitude of the distortion
# f1, f2, f3 : int, float, int
#     Filtering parameters.  f1 gives the size of the window, f2 gives the high-frequency cutoff, f3 gives the size of the filter
# xaxis_loc, yaxis_log : float
#     The locations to draw the x and y axes.  If not specified, they will be drawn from the bottom left of the plot
# xaxis_arrow, yaxis_arrow : str
#     Where to draw arrows on the x/y axes. Options are '+', '-', '+-', or ''
# ax_extend : float
#     How far (fractionally) to extend the drawn axes beyond the original axes limits
# expand_axes : bool
#     if True, then expand axes to fill the figure (useful if there is only a single axes in the figure)
def XKCDify(ax, mag=1.0,
            f1=50, f2=0.01, f3=15,
            bgcolor='w',
            xaxis_loc=None,
            yaxis_loc=None,
            xaxis_arrow='+',
            yaxis_arrow='+',
            ax_extend=0.1,
            expand_axes=False):
    # Get axes aspect
    ext = ax.get_window_extent().extents
    aspect = (ext[3] - ext[1]) / (ext[2] - ext[0])

    # Set X/Y axis position if not specified
    if xaxis_loc is None:
        xaxis_loc = ax.get_ylim()[0]

    if yaxis_loc is None:
        yaxis_loc = ax.get_xlim()[0]

    # Setup transformation
    trans = ax.transScale+ax.transLimits

    # Calculate x/y_min/max to have %5 free area around
    [[x_min, y_min]] = trans.inverted().transform(numpy.array([(-0.04*aspect, -0.04)]))
    [[x_max, y_max]] = trans.inverted().transform(numpy.array([(1+(0.04*aspect), 1.04)]))

    # Set the axis limits
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # # Draw axes
    [[yaxis_loc_in_axis_coordinates, xaxis_loc_in_axis_coordinates]] = trans.transform(numpy.array([(yaxis_loc, xaxis_loc)]))
    xaxis = pylab.Line2D([0.01, 0.99], [xaxis_loc_in_axis_coordinates, xaxis_loc_in_axis_coordinates], linestyle='-', color='k', solid_capstyle='round', transform=ax.transAxes)
    yaxis = pylab.Line2D([yaxis_loc_in_axis_coordinates, yaxis_loc_in_axis_coordinates], [0.01, 0.99] , linestyle='-', color='k', solid_capstyle='round', transform=ax.transAxes)

    # Label axes3, 0.5, 'hello', fontsize=14)
    ax.text(1.0-0.02, xaxis_loc_in_axis_coordinates-0.035, ax.get_xlabel(), fontsize=14, ha='right', va='top', rotation=6, transform=ax.transAxes)
    ax.text(yaxis_loc_in_axis_coordinates-0.03, 1-0.01, ax.get_ylabel(), fontsize=14, ha='right', va='top', rotation=78, transform=ax.transAxes)
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Add title
    ax.text(0.5, 1.0, ax.get_title(), ha='center', va='top', fontsize=18, transform=ax.transAxes)
    ax.set_title('')

    # Ticks
    ticks = []
    tick_size = 0.01
    tick_randomness = 0.0005

    for i in range(len(ax.get_xticks())):
        tick = ax.get_xticks()[i]
        text = ax.get_xticklabels()[i].get_text();

        if text == '':
            text = ax.xaxis.get_major_formatter().format_data_short(tick).strip()

        transform = ax.transAxes;

        [[tick_x, tick_y]] = trans.transform(numpy.array([(tick, xaxis_loc)]))

        x_int, y_int = xkcd_line(ax, [tick_x, tick_x], [tick_y+numpy.random.normal(tick_size, tick_randomness, 1)[0], tick_y-numpy.random.normal(tick_size, tick_randomness, 1)[0]], transform, mag, f1, f2, f3)
        line = pylab.Line2D(x_int, y_int, linestyle='-', color='k', solid_capstyle='round', transform=transform)
        ticks.append(line)       
        if (tick_x <= 1) and (tick_x >= 0):
            ax.text(tick_x, tick_y-numpy.random.normal(0.025, 0.001, 1)[0], text, fontsize=14, ha='center', va='top', rotation=-3+numpy.random.normal(1)*6, transform=ax.transAxes)

    for i in range(len(ax.get_yticks())):
        tick = ax.get_yticks()[i]
        text = ax.get_yticklabels()[i].get_text();

        if text == '':
            text = ax.yaxis.get_major_formatter().format_data_short(tick).strip()

        transform = ax.transAxes;

        [[tick_x, tick_y]] = trans.transform(numpy.array([(yaxis_loc, tick)]))

        x_int, y_int = xkcd_line(ax, [tick_x+numpy.random.normal(tick_size, tick_randomness, 1)[0], tick_x-numpy.random.normal(tick_size, tick_randomness, 1)[0]], [tick_y, tick_y], transform, mag, f1, f2, f3)
        line = pylab.Line2D(x_int, y_int, linestyle='-', color='k', solid_capstyle='round', transform=transform)
        ticks.append(line)
        if (tick_y <= 1) and (tick_y >= 0):
            ax.text(tick_x-numpy.random.normal(0.025, 0.001, 1)[0], tick_y, text, fontsize=14, ha='right', va='center', rotation=-3+numpy.random.normal(1)*6, transform=ax.transAxes)

    # Make all lines wiggly
    Nlines = len(ax.lines)
    lines = [xaxis, yaxis] + ticks + [ax.lines.pop(0) for i in range(Nlines)] 

    for line in lines:
        x, y = line.get_data()
        transform = line.get_transform()

        if line.get_linestyle() != 'None':
            x_int, y_int = xkcd_line(ax, x, y, transform, mag, f1, f2, f3)

            # create foreground and background line
            lw = line.get_linewidth()

            transform = line.get_transform();
            line.set_linewidth(2 * lw)
            line.set_data(x_int, y_int)

            # don't add background line for axes
            if (line is not xaxis) and (line is not yaxis) and (line not in ticks):
                line_bg = pylab.Line2D(x_int, y_int, color=bgcolor, linewidth=8 * lw, solid_capstyle='round', transform=transform)
                ax.add_line(line_bg)

        l1 = ax.add_line(line)

        # Something there a problems with more than 100 data points. Don't know why, but
        # http://matplotlib.1069221.n5.nabble.com/demo-curvelinear-grid-py-and-lines-over-100-points-td16153.html
        # suggests this as a workaround :)
        l1._transformed_path = None 
        l1._subslice = False 

    # Draw arrow-heads at the end of axes lines
    arr1 = 0.03 * numpy.array([-1, 0, -1])
    arr2 = 0.02 * numpy.array([-1, 0, 1])

    arr1[::2] += numpy.random.normal(0, 0.005, 2)
    arr2[::2] += numpy.random.normal(0, 0.005, 2)

    x, y = xaxis.get_data()
    if '+' in str(xaxis_arrow):
        ax.plot(x[-1] + arr1, y[-1] + arr2, color='k', lw=2, solid_capstyle='round', transform=ax.transAxes)
    if '-' in str(xaxis_arrow):
        ax.plot(x[0] - arr1, y[0] - arr2, color='k', lw=2, solid_capstyle='round', transform=ax.transAxes)

    x, y = yaxis.get_data()
    if '+' in str(yaxis_arrow):
        ax.plot(x[-1] + arr2, y[-1] + arr1, color='k', lw=2, solid_capstyle='round', transform=ax.transAxes)
    if '-' in str(yaxis_arrow):
        ax.plot(x[0] - arr2, y[0] - arr1, color='k', lw=2, solid_capstyle='round', transform=ax.transAxes)

    # Change all the fonts to humor-sans.
    prop = matplotlib.font_manager.FontProperties(fname='Humor-Sans.ttf', size=16)
    for text in ax.texts:
        text.set_fontproperties(prop)
    
    # modify legend
    leg = ax.get_legend()
    if leg is not None:
        leg.set_frame_on(False)
        
        for child in leg.get_children():
            if isinstance(child, pylab.Line2D):
                x, y = child.get_data()
                child.set_data(xkcd_line(ax, x, y, 'none', mag, f1, f2, f3));
                child.set_linewidth(2 * child.get_linewidth())
            if isinstance(child, pylab.Text):
                child.set_fontproperties(prop)
    
    if expand_axes:
        ax.figure.set_facecolor(bgcolor)
        ax.set_position([0, 0, 1, 1])
    
    ax.set_axis_off()

    return ax