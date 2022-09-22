install.packages("VennDiagram")

library("VennDiagram")

grid.newpage()

draw.triple.venn(area1 = 1186, area2 = 1876, area3 = 4120,
                 n12 = 616, n23 = 889, n13 = 765, n123 = 421,
                 category = c("E-cadherin", "HP1", "STAT92E"),
                 col = "Black", fill = c("Blue", "Magenta", "Green"))