# output 1
env1 = new.env()
load(here("output/stadem_results/LGR_STADEM_Chinook_2010.rda"),
     envir = env1)

# output 2
env2 = new.env()
output2 = load("C:/Git/SnakeBasinFishStatus/STADEM_results/LGR_STADEM_Chinook_2010.rda",
               envir = env2)

output1 = env1$stadem_mod
output2 = env2$stadem_mod

output1
output2
