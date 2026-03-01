# Counts Analysis


**Problem Statement**

When dealing with gene counts in an mRNA-seq dataset, it is important to normalize the data before performing any analyses so you can make accurate comparisons of gene expression between samples. There are a number of different methods you can use to accomplish and the method you use will depend on the kinds of samples you have and the analyses you want to perform. This first part of this assignment will guide you through two commonly used methods: Counts Per Million and DESeq2 Normalization.

Once you normalize your data you can do exploratory analysis and construct the appropriate graphs. As bioinformaticians, you will need to be able present your data in tidy reports with generated plots- sometimes between blocks of text. R Markdown is a good tool to accomplish this with, allowing you to display tables and plots alongside the bulk of your writing. You’ve had some exposure to R Markdown in previous assignments, now you will start writing your own code to generate the visualizations for your report. It will be important to remember that you need to keep your code separate depending on its function: the function declarations and implementation- your code that will work behind-the-scenes to perform data manipulation and process inputs without user interaction on the “back end”- should stay in your main.R file while your calls to construct and display the tibbles and visuals created in the back end should be put in your report.Rmd file- despite the user not directly manipulating or editing the outputs you generate, they are still interacting with the display by viewing it and therefore the lines of code that call the functions to create these outputs belong in your “front end.”


**Learning Objectives**
- Data normalization methods (CPM, DESeq2)
- Plotting data in Tidyverse
- Displaying results and compiling reports in R Markdown
- Application of Tidyverse functions


**Skill List**
- DESeq normalization, referencing the linked Bioconductor vignette if needed
- Tibble manipulation
- Creating different types of graphs in ggplot2
- Running PCA
