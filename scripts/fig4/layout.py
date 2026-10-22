"""Figure 4 layout: the size of every panel box and its position on the 183 mm page (all values in mm).

Every panel is drawn at exactly its box size by panels.py and assemble.py places the boxes 1:1, so the font
sizes set in panels.py are the printed sizes. Letters sit in a gutter left of their panel, except for the polar trees
(D, G), whose empty top-left corner holds the letter.

  row 1: A tree + character matrix | B spatial, clades
  row 2: C n_extant | D polar tree, fitness | E spatial, fitness
  row 3: F waterfall | G polar tree, fitness + Cldn4 | H spatial, Cldn4
"""

W = 183  # page width
M = 3  # page margin
LETTER = 4.5  # gutter for a panel letter
GAP = 2  # between panels in a row
ROW_GAP = 3.5  # between rows
H = 39.0  # every row

SPATIAL_ASPECT = 6481.2 / 3115.7  # width / height of the three tumor sections side by side (obsm["spatial"])
SPATIAL_W = 79.0  # width of the spatial strip in B, E, H (the colorbar is to its right)

# (letter, panel, width, height, letter gutter); rows are listed left to right
ROWS = [
    [("A", "a_tree_characters", 76.0, H, True), ("B", "b_spatial_clade", 90.0, H, True)],
    [("C", "c_n_extant", 36.0, H, True), ("D", "d_tree_fitness", 38.0, H, False), ("E", "e_spatial_fitness", 90.0, H, True)],
    [("F", "f_waterfall", 36.0, H, True), ("G", "g_tree_cldn4", 38.0, H, False), ("H", "h_spatial_cldn4", 90.0, H, True)],
]  # fmt: skip
BOX = {name: (w, h) for row in ROWS for _, name, w, h, _ in row}
# Placement offsets (dx, dy) in mm of a panel relative to its slot (positive dy = down); the letters stay in place
NUDGE = {"c_n_extant": (-1.5, 2.0), "g_tree_cldn4": (-1.5, 1.5)}
