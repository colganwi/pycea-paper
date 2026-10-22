"""Figure 3 layout: the size of every panel box and its position on the 183 mm page (all values in mm).

Every panel except the polar tree (A) is drawn at exactly its box size by panels_data.py / panels_schematics.py, and
assemble.py places the boxes 1:1, so the font sizes set in the panel scripts are the printed sizes. Panel letters sit
on the top-left corner of their box (the panel scripts keep that corner empty), except G, which has its own gutter.

  row 1: A polar tree | B UMAP (germ layer) | C UMAP (lineage)
  row 2: D schematic + UMAP of neural crest linkage | E schematic + pairwise linkage heatmap
  row 3: F Moran's I schematic + ranked heritability | G UMAPs of the four most heritable genes
"""

W = 183  # page width
M = 3  # page margin
LETTER = 4.5  # gutter for a letter that does not sit on its box (G)
SUB_GAP = 1.5  # between the parts of one panel (schematic | data plot)
GAP = 2  # between panels in a row
ROW_GAP = 3.5  # between rows

# Shared drawing sizes
UMAP_ASPECT = 38.7 / 43.7  # width / height of the UMAP data range (data/fig3_cells.csv.gz)
UMAP_W = 46.0  # width of the B and C UMAPs
TREE_W, TREE_H = 28.0, 18.0  # every tree schematic (D, E, F)
TITLE_X = 4.8  # x of schematic titles, right of the panel letter

# Row heights
H1, H2, H3 = 60.5, 62.0, 59.0

# (letter, panel, width, height, letter gutter); rows are listed left to right
ROWS = [
    [("A", "a_tree", 53.5, None, False), ("B", "b_umap_germ_layer", 46.5, H1, False), ("C", "c_umap_lineage", 73.0, H1, False)],
    [("D", "d_schematic", 41.5, H2, False), (None, "d_umap", 42.5, H2, False), ("E", "e_schematic", 40.0, H2, False),
     (None, "e_heatmap", 48.0, H2, False)],
    [("F", "f_schematic", 44.0, H3, False), (None, "f_waterfall", 62.5, H3, False), ("G", "g_umaps", 62.5, H3, True)],
]  # fmt: skip
BOX = {name: (w, h) for row in ROWS for _, name, w, h, _ in row}
