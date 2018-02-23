options(finemapr_caviar = "~/apps/paintor/CAVIAR")

ex <- example_finemap()

out <- run_paintor(ex$tab1, ex$ld1, ex$n1, args = "-enumerate 3")

print(out)

plot(out)
