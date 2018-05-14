#options(finemapr_caviar = "~/apps/paintor/CAVIAR")

ex <- example_finemap()

#out <- run_paintor(ex$tab1, ex$ld1, ex$n1, args = "-enumerate 3")
out <- finemapr(list(ex$tab1, ex$tab2), list(ex$ld1, ex$ld2), list(ex$n1, ex$n2), method = "paintor", args = "-enumerate 4")

print(out)

plot(out)
