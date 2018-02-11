
ex <- example_finemap()

out <- run_finemap(ex$tab1, ex$ld1, ex$n1, args = "--n-causal-max 4")

print(out)

plot(out)
