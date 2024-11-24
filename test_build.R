library(styler)
style_dir(path = "/Users/zhiz/zz/github/DirichletMultinomial", filetype = c("R", "Rmd"))

# rmarkdown::render("README.Rmd", output_format = "md_document")

x <- rmarkdown::render("README.Rmd", run_pandoc = FALSE, clean = FALSE)
knit_meta <- attr(x, "knit_meta")
rmarkdown::render(input = "README.md", knit_meta = knit_meta)
