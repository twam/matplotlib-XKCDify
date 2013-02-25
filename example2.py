# Example from http://jakevdp.github.com/blog/2012/10/07/xkcd-style-plots-in-matplotlib/ in Action

import numpy
import pylab
import XKCDify

# Some helper functions
def norm(x, x0, sigma):
    return numpy.exp(-0.5 * (x - x0) ** 2 / sigma ** 2)

def sigmoid(x, x0, alpha):
    return 1. / (1. + numpy.exp(- (x - x0) / alpha))
    
# define the curves
x = numpy.linspace(0, 1, 100)
y1 = numpy.sqrt(norm(x, 0.7, 0.05)) + 0.2 * (1.5 - sigmoid(x, 0.8, 0.05))
y2 = 0.2 * norm(x, 0.5, 0.2) + numpy.sqrt(norm(x, 0.6, 0.05)) + 0.1 * (1 - sigmoid(x, 0.75, 0.05))
y3 = 0.05 + 1.4 * norm(x, 0.85, 0.08)
y3[x > 0.85] = 0.05 + 1.4 * norm(x[x > 0.85], 0.85, 0.3)

# draw the curves
ax = pylab.axes()
ax.plot(x, y1, c='gray')
ax.plot(x, y2, c='blue')
ax.plot(x, y3, c='red')

# Add some text
ax.text(0.3, -0.1, "Yard")
ax.text(0.5, -0.1, "Steps")
ax.text(0.7, -0.1, "Door")
ax.text(0.9, -0.1, "Inside")

# Add some annotations
ax.text(0.05, 1.1, "fear that\nthere's\nsomething\nbehind me")
ax.plot([0.15, 0.2], [1.0, 0.2], '-k', lw=0.5)

ax.text(0.25, 0.8, "forward\nspeed")
ax.plot([0.32, 0.35], [0.75, 0.35], '-k', lw=0.5)

ax.text(0.8, 0.6, "embarrassment")
ax.plot([0.95, 0.8], [0.65, 1.05], '-k', lw=0.5)

# Set title
ax.set_title("Walking back to my\nfront door at night:")

# Set limits
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.5)

# No x/yticks
ax.set_xticks([])
ax.set_yticks([])

# modify all the axes elements in-place
XKCDify.XKCDify(ax, expand_axes=True)

# save to file
pylab.savefig('example2.pdf')
