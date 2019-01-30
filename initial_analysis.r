initial <- read.table("initial.out")

print("Instances:")
print(nrow(initial)) #Total number of instances
print("Initial analysis wins vs grevlex")
print(nrow(initial[initial$V10 <= initial$V4, ])) #When initial analysis beats grevlex
print("Initial analysis wins vs few iterations of Perry")
print(nrow(initial[initial$V10 <= initial$V7, ])) #When initial analysis beats dynamic, stopped at m iterations
