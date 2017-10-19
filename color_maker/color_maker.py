#!/usr/bin/env python


def color_map(name=None):
    """ Function to load any of the NCL colormaps as
    matplotlib LinearSegmented Colormaps.  A list of all the color
    tables available is at:
        http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml
    Most color tables are from NCL.

    To use:
        cmap = color_map('NAME OF COLORMAP HERE')
        e.g. cmap = color_map('lithology')
    """

    cmap_path = '/home/disk/meso-home/bzim/relampago/color_maker/colormaps'
    from matplotlib.colors import LinearSegmentedColormap
    import os
    if name + '.coltbl' in os.listdir(cmap_path):
        infile = open('/'.join((cmap_path, name+'.coltbl')),'r')
    elif name + '.rgb' in os.listdir(cmap_path):
        infile = open('/'.join((cmap_path, name+'.rgb')),'r')
    else:
        print("Unable to find colormap:", name)
        return None

    # Build the colors
    lines = infile.readlines()
    use_lines = []
    maxval = 0
    for line in lines:
        if line[0].isalpha() or line[0] == '#':
            continue
        linesp = line.split()
        if len(linesp) < 3:
            continue
        fline = [float(x) for x in linesp]
        maxval = max(fline + [maxval])
        use_lines.append(fline)
    lines = use_lines
    nlines = len(lines)
    if maxval > 1.0:
        div = 255.
    else:
        div = 1.

    infile.close()
    cdict = {}
    cdict['red'] = [(n/float(nlines-1), line[0]/div, line[0]/div) for n,line in enumerate(lines)]
    cdict['green'] = [(n/float(nlines-1), line[1]/div, line[1]/div) for n,line in enumerate(lines)]
    cdict['blue'] = [(n/float(nlines-1), line[2]/div, line[2]/div) for n,line in enumerate(lines)]

    return LinearSegmentedColormap(name=name, segmentdata=cdict)


if __name__ == '__main__':
    color_map('radar_1')
