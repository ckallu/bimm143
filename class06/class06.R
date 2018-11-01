#' ---
#' title: "Class 6"
#' output: github_document
#' ---

# Reading input data as three separte tables

read.csv("test1.txt", header = TRUE)

read.table("test2.txt", header = TRUE, sep = "$")

read.table("test3.txt")

# 1A

# Code with duplication and errors

df <- data.frame(a=1:10, b=seq(200,400,length=10),c=11:20,d=NA)
df$a <- (df$a - min(df$a)) / (max(df$a) - min(df$a))
df$b <- (df$b - min(df$a)) / (max(df$b) - min(df$b))
df$c <- (df$c - min(df$c)) / (max(df$c) - min(df$c))
df$d <- (df$d - min(df$d)) / (max(df$a) - min(df$d)) 

df

# function to scale

scaling <- function (x) {
  rng = range (x)
  x = (x - rng[1] / rng[2] - rng[1])
  return(x)
}

# use function to scale columns in df 

df$a <- scaling (df$a)

df$b <- scaling (df$b)

df$c <- scaling (df$c)

df$d <- scaling (df$d)

df

#1B

library("bio3d")
# Function that takes the name of a pdb file as input and analyzes protein drug interactions for that protein and outputs a plot for the specified protein.

protein_analysis <- function (x) {
  
  # s stores the info obtained from reading the pdb file.
  s <- read.pdb(x)
  
  s.chainA <- trim.pdb(s, chain = "A", elety = "CA")
  s.b <- s.chainA$atom$b
  plot_protein <- plotb3(s.b, sse = s.chainA, typ = "l", ylab = "Bfactor", main = x)
  
  #outputs a plot showing protein drud interaction for the specified protein
  return(plot_protein)
}

protein_analysis("4AKE")
protein_analysis("1AKE")
protein_analysis("1E4Y")

