# Example from http://jakevdp.github.com/blog/2012/10/07/xkcd-style-plots-in-matplotlib/ in Action

import numpy
import pylab
import XKCDify

ax = pylab.axes()

x = numpy.linspace(0, 10, 100)
ax.plot(x, numpy.sin(x) * numpy.exp(-0.1 * (x - 5) ** 2), 'b', lw=1, label='damped sine')
ax.plot(x, -numpy.cos(x) * numpy.exp(-0.1 * (x - 5) ** 2), 'r', lw=1, label='damped cosine')

# set labels/title
ax.set_title('check it out!')
ax.set_xlabel('x label')
ax.set_ylabel('y label')

# set legend position
ax.legend(loc='lower right')

# set x/y limits
ax.set_xlim(0, 10)
ax.set_ylim(-1.0, 1.0)

# No x/yticks
ax.set_xticks([])
ax.set_yticks([])

# modify all the axes elements in-place
XKCDify.XKCDify(ax, xaxis_loc=0.0, yaxis_loc=1.0, xaxis_arrow='+-', yaxis_arrow='+-', expand_axes=True)

# save to file
pylab.savefig('example.pdf')
