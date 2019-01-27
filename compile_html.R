library(rmarkdown)

rmarkdown::render(
		  "Rmd/GWAS_practicals.Part_1.Rmd", 
		  "html_document", 
		  output_dir='html')

rmarkdown::render(
                  "Rmd/GWAS_practicals.Part_2.Rmd",
                  "html_document",
                  output_dir='html')

rmarkdown::render(
                  "Rmd/GWAS_practicals.Part_3.Rmd",
                  "html_document",
                  output_dir='html')
