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
