source("http://tinyurl.com/rescale-R")

rescale2 ( c(1,10,"string"))

# Lets define an example x and y
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)

is.na(x) & is.na(y)

sum (is.na(x))
sum (is.na(y))

sum ((is.na(x) & is.na(y)))

# function

both_na <- function (x,y) {
  sum ((is.na(x) & is.na(y)))
}

x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA, NA, NA)

both_na(x,y1)
both_na(x,y2) # not what we want here!

both_na2 <- function (x,y) {
  if(length(x) != length(y)){
    stop("Input x and input y should be of the same length")
  }
  sum ((is.na(x) & is.na(y)))
}

both_na2(x,y2)

both_na3(x,y1)

# source("http://tinyurl.com/rescale-R")
# Start with a simple version of the problem
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
                  exp=c(2,1,1),
                  stringsAsFactors=FALSE)
df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
                  exp=c(-2, NA, 1, 2),
                  stringsAsFactors=FALSE)
# Simplify further to single vectors
x <- df1$IDs
y <- df2$IDs

intersect(x,y)
x %in% y

gene_intersect(x,y)


source("https://bioconductor.org/biocLite.R")
install.packages('blogdown')
library(blogdown)
install_hugo()
blogdown::new_site()
