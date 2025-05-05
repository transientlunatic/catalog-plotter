import matplotlib.colors
from cycler import cycler

colours = """#3D49E3
#74E33D
#E36F3D
#8E6756
#88C6D4
#C88EC7
#48070E
#FC851D
#4DCA70
#5BB75A
#0D2421
#A91031
#3062D3
#D7E83B
#741D0A
#BE8FCD""".split()


colours = [matplotlib.colors.to_rgba(color) for color in colours]
    
color_cycle = cycler('color', list(colours))

