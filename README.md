# Reproducible research: version control and R

### Question 4

**a.**  
The .R script provided contains the code to write the random_walk() function. The random_walk() function creates a data frame and fills it with positions based on randomly generated angles for multiple iterations. This simulates the random walk. random_walk was run twice and produces two data frames (data1 and data2), each with 500 steps (i.e. 500 iterations). 

The walks are visualised graphically by plotting both data1 and data2. The grid.arrange() function places the two graphical visualisations of the two walks next to each other. The graphs show the walks in two dimensions (in a grid with a x and y-axis). The direction of the path is indicated by the colour gradient; dark blue represents when t is small/the earlier steps, and light blue represents when t is larger/the later steps. Therefore, the path moves from dark blue to light blue.

The random_walk() function produces a new data frame each time, i.e. re-running the random_walk(500) makes a different data1 and data2 and therefore different paths each time. Executing the whole code repeatedly creates new, random pairs of graphs each time.

**b.**  
A random seed is the starting point for generating random numbers in functions in R. When R generates random numbers, the numbers are not truly random but instead are pseudorandom numbers. It is pseudorandom because the algorithm that is used to generate the randomness uses a ‘seed’ to initialise it. The seed is randomly selected and therefore a random outcome is produced. This is stored in .Random.seed in the global environment and is a vector of integers (Random-R documentation, n.d). 

However, knowing this seed means that the random numbers can be predicted and therefore reproduced. Setting the same random seed means that the same sequence of random numbers is generated each time. Crucially, setting seed is important for reproducibility because the exact same results–despite being random–are produced each time. Therefore, the code can be run by another person and they can achieve the exact same outcome. Setting a random seed is important for other functions, e.g. debugging a code (setting a seed would make it easier to identify errors by removing the element of randomness) (r-coder.com, n.d).

The random seed can be set with the function: set.seed(n) from the base R package. Sharing the value, n, used in set.seed(n) allows random number generation to be reproduced.

**c.**  
Please see the script and part d. for the edits. 

**d.**  
Below are images of the latest commit in the comparison view.  
<p align="center">
   <img width="1434" src="https://github.com/pepperepperepper/reproducible-research_homework/blob/main/images/q4_partd_image1.png">
</p>
<p align="center">
   <img width="1434" src="https://github.com/pepperepperepper/reproducible-research_homework/blob/main/images/q4_partd_image2.png">
</p>

***

### Question 5  

**a.**  
The dsDNA virus data has 33 rows and 13 columns.

The column names (before cleaning and converting to snake case) are: “Family”, “Genus”, “Type.species”, “GenBank.accession.no.”, “Envelope”, “Virion.type”, “T”, “Virion.diameter..nm.”, “Virion.length..nm.”, “Virion.volume..nm.nm.nm.”, “Molecule”, “Genome.length..kb.”, and “Protein.no.”

**b.**  
A log transformation can be performed on both “genome_length_kb” and “virion_volume_nm_nm_nm”. A log transformation is applicable here because both genome lengths and virion volume span several magnitudes (i.e. is exponential). Converting to a log scale linearises the relationship. As a result, the values for genome length range from roughly 1.5 to 8 and values for virion volume ranges from roughly 10 to 18. The linear relationship assumption of linear regression is now met and a linear model can be fitted.

The log() function carries out the log transformation and is applied to the whole column of data (dsvirus$genome_length_kb and dsvirus$virion_volume_nm_nm_nm). 

**c.**  
After linearising the relationship by taking logs, a linear model can be fitted using the lm() function. Using summary() prints the results of the model and coef() shows the y-intercept and the slope of the line.  

The results of the model:
```math
\begin{equation}
y-intercept = 7.074800
\end{equation}
```
This has a p-value of 2.28e-10.

```math
\begin{equation}
slope = 7.074800
\end{equation}
```
This has a p-value of 6.44e-10.

The p-values for the y-intercept and slope are both vanishingly small and, at a 0.05 significance level, we can therefore conclude that these results are statistically significant.

The equation of the line is therefore:
```math
\begin{equation}
V = 7.074800 + 1.515228L
\end{equation}
```
Where V = virion volume and L = length of genome.

But we took logs of both genome length and virion volume, so the equation is:

```math
\begin{equation}
log(V) = 7.074800 + 1.515228*log(L)
\end{equation}
```

Both sides of the equation need to be exponentiated to reverse the logs:
```math
\begin{equation}
e^{log(V)} = e^{7.074800 + 1.515228*log(L)}
\end{equation}
```

This simplifies to:
```math
\begin{equation}
V = e^{7.074800 + 1.515228*log(L)}
\end{equation}
```

Due to the power rule, this can be rearranged to:
```math
\begin{equation}
V = e^{7.074800} + e^{1.515228*log(L)}
\end{equation}
```

And finally simplified again:
```math
\begin{equation}
V = e^{7.074800} + L^{1.515228}
\end{equation}
```

The y-intercept can be calculated:
```math
\begin{equation}
V = 1181.807 * L^{1.515228}
\end{equation}
```

Therefore,

```math
\begin{equation}
β = 1181.807
\end{equation}
```
```math
\begin{equation}
α = 1.515228
\end{equation}
```

Where α = allometric exponent, and β = scaling factor.  

The values found in table 2 for dsDNA viruses were:
```math
\begin{equation}
α = 1.52
\end{equation}
```
```math
\begin{equation}
β = 1182
\end{equation}
```
This is consistent with the results of the data analysis; the results in the paper are rounded to 2 d.p for α and to the nearest integer for β.

**d.**  
Please find below the code that reproduces the figure below (from setting up the workspace -> cleaning -> transforming -> plotting).  
```{r reproduce-figure-code}
# Install and load packages
install.packages(c("janitor", "ggplot2", "dplyr", "magrittr"))
library(janitor)
library(ggplot2) 
library(dplyr)
library(magrittr)

# Load in the data
dsvirus <- read.csv("Cui_etal2014.csv")

# View the data
names(dsvirus) # view names of the columns
head(dsvirus) # view first few rows of the data

# Clean the data
clean_dsvirus <- dsvirus %>%
  clean_names() %>% # convert to lower case and snake case
  select(c("genome_length_kb", "virion_volume_nm_nm_nm")) # select only the relevant columns

# Log transform the data to linearise the relationship
clean_dsvirus$log_genome_length <- log(clean_dsvirus$genome_length_kb)
clean_dsvirus$log_virion_volume <- log(clean_dsvirus$virion_volume_nm_nm_nm)

# Plot the data
plot_dsvirus <- ggplot(clean_dsvirus, aes(x = log_genome_length, y = log_virion_volume)) +
    geom_point(color = "black",
               size = 2,
               alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    labs(x = "log[Genome length(kb)])",
         y = "log[Virion volume (nm³)]",) +
    theme_bw()
print(plot_dsvirus)

```

**e.**
The equation is:

```math
\begin{equation}
V = 1181.807 * L^{1.515228}
\end{equation}
```

For a 300 kb dsDNA virus, L = 300.  
Therefore, by we can sub L = 300 into the equation:  

```math
\begin{equation}
V = 1181.807 * 300^{1.515228}
\end{equation}
```

```math
\begin{equation}
V = 6698076
\end{equation}
```  

Therefore, the estimated volume, V for a 300 kb dsDNA virus is **6,698,076nm³**.

**References**

r-coder.com, 'Setting the Seed in R for Reproducibility'. 
Available at: https://r-coder.com/set-seed-r/ (Accessed: 6 December 2023).

R Project for Statistical Computing, 'Random - R Documentation'. Available at: https://stat.ethz.ch/R-manual/R-devel/library/base/html/Random.html (Accessed: 6 December 2023).


## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points (plus an optional bonus question worth 10 extra points). First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers. All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   - A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)
   - Investigate the term **random seeds**. What is a random seed and how does it work? (5 points)
   - Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points)
   - Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points)

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \beta L^{\alpha}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   - Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)
   - What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)
   - Find the exponent ($\alpha$) and scaling factor ($\beta$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)
   - Write the code to reproduce the figure shown below. (10 points)

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  - What is the estimated volume of a 300 kb dsDNA virus? (4 points)

**Bonus** (**10 points**) Explain the difference between reproducibility and replicability in scientific research. How can git and GitHub be used to enhance the reproducibility and replicability of your work? what limitations do they have? (e.g. check the platform [protocols.io](https://www.protocols.io/)).
