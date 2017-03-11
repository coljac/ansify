#!env python
# -*- coding: utf-8 -*-
# TODO: Aspect ratio messedup
# TODO: Color bug - see trump approx vs full
# TODO: The white "fringes"
# TODO: Three-quarters and other figures
"""
Huh?
"""
import sys
import numpy as np
import os
import colors
from PIL import Image
import coltools as ct
from colors import color

import colormath
from colormath.color_objects import LabColor
from colormath.color_diff import delta_e_cie1976
from colormath.color_objects import XYZColor, sRGBColor
from colormath.color_conversions import convert_color
from rgbtoshort import rgb2short

# ░▒▓
#▐░▒▓▔▕▖▗▘▙▚▛▜▝▞▟

FULL_BLOCK = "█"
RIGHT_HALF_BLOCK = "▐"
LEFT_HALF_BLOCK = "▌"
UPPER_HALF_BLOCK = "▀"
LOWER_HALF_BLOCK = "▄"

colormap = None
diffs = {}

# Possibilities for color:
# Full block
# Quadrant 1, 2, 3, 4
# Top half, bottom half, right half, left half
# 1 and 4, 2 and 3
color_pairs = ct.fat("colors.txt")
ansicolors = np.load('ansicolors.npz')['colors']

def color_difference(color1, color2):
    c1 = convert_color(sRGBColor(*color1), LabColor, target_illuminant='d50')
    c2 = convert_color(sRGBColor(*color2), LabColor, target_illuminant='d50')
    # c1 = LabColor(lab_l=0.9, lab_a=16.3, lab_b=-2.22)
    # c2 = LabColor(lab_l=0.7, lab_a=14.2, lab_b=-1.80)
    delta_e = delta_e_cie1976(c1, c2)
    return delta_e


def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def get_best_ansi_color(rgb, mode="full"):
    if mode == "full":
        return get_best_ansi_color_full(rgb)
    elif mode == "approx":
        return get_best_ansi_color_approx(rgb)
    elif mode == "quick":
        return get_best_ansi_color_quick(rgb)
    elif mode == "rgb":
        return get_best_ansi_color_rgb(rgb)

def get_best_ansi_color_approx(rgb):
    global colormap
    if colormap is None:
        print("Loading color map.")
        colormap = np.load('colormap_approx.npz')['colormap'].astype(int)
    a, b, c = [int(min(24, rgb[a]/10)) for a in range(3)]
    approx = colormap[a, b, c]
    return int(approx)

def get_best_ansi_color_quick(rgb):
    hexv = rgb_to_hex(*rgb)
    vals = rgb2short(hexv)
    return int(vals[0])

def get_best_ansi_color_rgb(rgb):
    rgb_color_dist = color_distances(rgb)
    closest = np.where(rgb_color_dist == rgb_color_dist.min())[0]
    return int(closest[0])

def color_distances(rgb):
    global ansicolors
    distances = np.zeros(256)
    distances = np.power(ansicolors[:, 0] - rgb[0], 2) + np.power(ansicolors[:, 1] - rgb[1], 2) + \
        np.power(ansicolors[:, 2] - rgb[2], 2)
    distances = np.sqrt(distances)
    return distances

def get_best_ansi_color_full(rgb):
    global colormap
    global color_pairs
    global ansicolors

    if colormap is None:
        print("Loading color map.")
        colormap = np.load('colormap_full.npz')['colormap'].astype(int)

    rgb = rgb.astype(np.uint16)
    best_match = colormap[rgb[0], rgb[1], rgb[2]]
    if best_match < 999:# and False:
        return int(best_match)

    lowest_distance = 1.0e9
    best_match = 0
    max_distance = 70

    rgb_color_dist = color_distances(rgb)
    in_range = np.where(rgb_color_dist < max_distance)[0]

    for j in range(len(in_range)):
        i = in_range[j]
        ansi_color = ansicolors[i]
        distance = color_difference(rgb, ansi_color)
        if distance < lowest_distance:
            lowest_distance = distance
            best_match = i
    rgb = rgb.astype(np.uint16)
    colormap[rgb[0], rgb[1], rgb[2]] = best_match
    return int(best_match)

def main():
    """ I need to fix the linter """

    infile = sys.argv[1]
    outfile = infile.split(".")[0] + ".ans"

    if len(sys.argv) > 2:
        outfile = sys.argv[2]
        
    width = 60
    if len(sys.argv) > 3:
        width = int(sys.argv[3])
    mode = "full"
    if len(sys.argv) > 4:
        mode = sys.argv[4]
    height = int(width * 0.47)
    size = height, width

    with open(outfile, "w") as out:
        try:
            # Assume for now square images. TODO

            im = Image.open(infile)
            # im.thumbnail(size, Image.ANTIALIAS)
            im = im.resize((size[1]*2, size[0]*2), Image.ANTIALIAS) # wxh
            # im.save("thumbnail.png")
            data = np.asarray(im, dtype=np.uint16)
            # print("Resized pixel shape:", data.shape)
            done = 0
            total = size[0] * size[1]
            for y in range(size[0]): # height
                for x in range(size[1]): # width
                    pixels = data[y*2:(y*2)+2, x*2:(x*2)+2]
                    distance = 0
                    top = pixels[0, :]
                    bottom = pixels[1, :]
                    left  = pixels[:, 0]
                    right = pixels[:, 1]
                    right_diag = ((pixels[0, 1] + pixels[1, 0])/2)#.astype(int)
                    left_diag = ((pixels[0, 0] + pixels[1, 1])/2)#.astype(int)
                     
                    dist_topbot = color_difference(top.mean(axis=0), bottom.mean(axis=0))
                    dist_lr = color_difference(left.mean(axis=0), right.mean(axis=0))
                    dist_diag = color_difference(left_diag, right_diag)

                    dist_1 = color_difference(pixels[0, 0], (pixels.sum(axis=0).sum(axis=0) - pixels[0, 0])/3)
                    dist_2 = color_difference(pixels[1, 0], (pixels.sum(axis=0).sum(axis=0) - pixels[1, 0])/3)
                    dist_3 = color_difference(pixels[0, 1], (pixels.sum(axis=0).sum(axis=0) - pixels[0, 1])/3)
                    dist_4 = color_difference(pixels[1, 1], (pixels.sum(axis=0).sum(axis=0) - pixels[1, 1])/3)

                    distances = np.array([dist_lr, dist_topbot, dist_diag, dist_1, dist_2,
                        dist_3, dist_4])
                    contrast = np.argmax(distances)

                    if contrast == 0:
                        best_color_left = get_best_ansi_color(left.mean(axis=0), mode=mode)
                        best_color_right = get_best_ansi_color(right.mean(axis=0), mode=mode)
                        if best_color_right == best_color_left:
                            out.write(color(" ", bg=best_color_left))
                        else:
                            out.write(color("▐", bg=best_color_left, fg=best_color_right))
                    elif contrast == 2:
                        best_color_left = get_best_ansi_color(left_diag, mode=mode)
                        best_color_right = get_best_ansi_color(right_diag, mode=mode)
                        if best_color_right == best_color_left:
                            out.write(color(" ", bg=best_color_left))
                        else:
                            out.write(color("▞", bg=best_color_left, fg=best_color_right))
                    elif contrast == 1:
                        best_color_top = get_best_ansi_color(top.mean(axis=0), mode=mode)
                        best_color_bottom = get_best_ansi_color(bottom.mean(axis=0), mode=mode)
                        if best_color_bottom == best_color_top:
                            out.write(color(" ", bg=best_color_top))
                        else:
                            out.write(color("▄", bg=best_color_top, fg=best_color_bottom))
                    else:
                    #▐░▒▓▔▕▖▗▘▙▚▛▜▝▞▟
                    # dist_1 = color_difference(pixels[0, 0], (pixels.sum(axis=0).sum(axis=0) - pixels[0, 0])/3)
                    # dist_2 = color_difference(pixels[1, 0], (pixels.sum(axis=0).sum(axis=0) - pixels[1, 0])/3)
                    # dist_3 = color_difference(pixels[0, 1], (pixels.sum(axis=0).sum(axis=0) - pixels[0, 1])/3)
                    # dist_4 = color_difference(pixels[1, 1], (pixels.sum(axis=0).sum(axis=0) - pixels[1, 1])/3)
                        shapes = ("▟","▜","▙","▛")
                        shape = shapes[contrast-3]
                        if contrast == 3:
                            single = pixels[0, 0]
                        elif contrast == 4:
                            single = pixels[1, 0]
                        elif contrast == 5:
                            single = pixels[0, 1]
                        else:
                            single = pixels[1, 1]
                        rest = (pixels.sum(axis=0).sum(axis=0) - single)/3

                        best_color_single = get_best_ansi_color(single, mode=mode)
                        best_color_rest = get_best_ansi_color(rest, mode=mode)
                        if best_color_single== best_color_rest:
                            out.write(color(" ", bg=best_color_single))
                        else:
                            out.write(color(shape, bg=best_color_single, fg=best_color_rest))
                    done += 1
                    ct.progbar(done, total)
                out.write("\n")

            if mode == "full":
                np.savez("colormap_full.npz", colormap=colormap)

        except IOError, e:
            print("cannot create thumbnail for '%s'" % infile)
 

if __name__ == "__main__":
    main()

