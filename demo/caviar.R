options(finemapr_caviar = "~/apps/caviar/CAVIAR")

ex <- example_finemap()

out <- run_caviar(ex$tab1, ex$ld1, args = "-c 3")

print(out)

plot(out)
