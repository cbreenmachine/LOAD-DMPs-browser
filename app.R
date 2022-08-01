
library(BiocManager)
options(repos = BiocManager::repositories())
library(data.table)
library(tidyverse)
library(ggbio)
library(EnsDb.Hsapiens.v86)
library(cowplot)
library(ggsci)
library(shiny)
library(scales)


global_size = 24

# show_col(pal_jama("default", alpha=0.95)(7))
colors <- pal_jama("default", alpha=1)(4)
# show_col(colors)

c_gene <- colors[4]
c_sig <- colors[2]
c_not_sig <- colors[3]


#--> Load in data (DO NOT CHANGE NAMES)
DMPs.df <- fread("app.DMPs.csv") %>% dplyr::mutate(chrom = factor(chrom))
genes.df <- fread("app.genes.csv") %>% dplyr::mutate(gene = factor(gene), chrom = factor(chrom))
ensdb <- EnsDb.Hsapiens.v86 # DO NOT CHANGE



# Theme -------------------------------------------------------------------
theme_dmp = function(){
  theme_minimal() %+replace%
  theme(strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.9,.87),
        text = element_text(size=global_size),
        plot.title = element_text(hjust = 0.5)
  ) 
}


# Gene model function -----------------------------------------------------
plot_gene_model <- function(ichrom, igene, strand='*', db=ensdb){
  # Top of plot, plots the gene introns/exons
  iseq <- str_remove(ichrom, pattern='chr0|chr')
  my_filter <- AnnotationFilter(~ symbol == igene & seq_name == iseq)
  
  p <- autoplot(db, my_filter, 
                stat="reduce", label=F, 
                color= c_gene, fill= c_gene) +
    ggtitle(paste0(igene, " gene model (", strand, " strand)")) +
    ylab("") +
    scale_x_continuous(position = "top") +
    theme_dmp()
    
  return(p@ggplot)
}

# Not, is
sig_levels <- c("LFDR \u2265 0.05", "LFDR < 0.05")
#--> Bottom of plot, scatter of 
plot_points <- function(chrom, start, end, df=DMPs.df){
  p <- df %>%
    dplyr::filter(chrom == chrom) %>%
    dplyr::filter(position > start & position + 1 < end) %>%
    dplyr::mutate(Significance = 
                    factor(
                      ifelse(abs(signed.logged.lfdr) <  -log10(0.05),
                           sig_levels[1], sig_levels[2]),
                      levels = sig_levels)) %>%
    ggplot(aes(x = position, y = signed.logged.lfdr, color = Significance)) +
    geom_point(size = 3) +
    xlab("Genomic position") +
    scale_color_manual(values = c(c_sig, c_not_sig)) +
    xlim(c(start, end)) +
    ylim(c(-3, 3)) +
    geom_hline(yintercept = 0, color = "black", alpha = 0.5) +
    geom_hline(yintercept = log10(0.05), color = "black", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), color = "black", alpha = 0.3) +
    ylab(expression('-log'[10]*'(LFDR)')) +
    theme_dmp()
  return(p)
}



join_via_cow <- function(p1, p2){
  p <- cowplot::plot_grid(p1, p2, ncol = 1, axis = "b",
                          align = "v", rel_heights = c(1, 2)) +
    theme(plot.background = element_rect(fill = "white", colour = "white"))
  return(p)
}



#--> Driver function
make_figure <- function(ichrom, igene, left=NULL,right=NULL, df=genes.df, save=FALSE){
  row <- df %>% dplyr::filter(chrom==ichrom, gene == igene)
  print(row)
  chrom <- dplyr::pull(row, 'chrom')
  start <- dplyr::pull(row, 'start')
  end <- dplyr::pull(row, 'end')
  strand <- dplyr::pull(row, 'strand')
  
  #TODO: add left/right sliders
  p1 <- plot_gene_model(ichrom, igene, strand=strand)
  p2 <- plot_points(chrom, start, end)
  
  if (!is.null(left) & !is.null(right)){
    print("Custom x lims")
    p1 <- p1 + xlim(c(left, right))
    p2 <- p2 + xlim(c(left, right))
  }
  p <- join_via_cow(p1, p2)
  
  if (save){
    cowplot::save_plot(filename=paste0(igene, ".png"), plot = p, 
                       base_width = 9, base_height = 6)
  }
  
  
  return(p)
}



start_gene <- 'BRCA1'
row <- genes.df %>% dplyr::filter(gene == start_gene)
start_left <- dplyr::pull(row, start)
start_right <- dplyr::pull(row, end)
start_chrom <- dplyr::pull(row, chrom)


ui <- pageWithSidebar(
  headerPanel("Differentially Methylated Positions"),
  sidebarPanel(
    selectizeInput(inputId ='chrom', 'Select chromosome', choices = c("choose" = "", levels(genes.df$chrom)), selected=start_chrom),
    selectizeInput(inputId = 'gene', 'Select gene', choices = c("choose" = "", levels(genes.df$gene)), selected=start_gene),
    sliderInput(inputId = 'position', min=start_left, max=start_right, value = c(start_left, start_right), label = "Genomic position"),
    actionButton(inputId = 'update', "Generate Plot"),
    
    shiny::br(),

    shiny::p("Plots take a few seconds to render of clicking 'Generate Plot'. 
             Data shown here include 42 LOAD versus 42 controls."),
    shiny::br(),
    shiny::p("For feature requests and bug reports, please email cebreen@wisc.edu")
  ),
  
  mainPanel(
    plotOutput(outputId = "geneModelPlot", width = 1200, height = 900)
  )
)




server <- function(input, output, session) {
  # Handle filtering "externally"
  observeEvent(input$chrom, {
    updateSelectizeInput(
      session, 
      "gene", 
      choices = genes.df %>% dplyr::filter(chrom == input$chrom) %>% dplyr::pull(gene) %>% sort()
    )
  }, ignoreInit = T)


  
  observeEvent(input$gene, {
    row <- genes.df %>% dplyr::filter(chrom == input$chrom, gene == input$gene)
    updateSliderInput(session,
                      inputId = "position",
                      min=row$start,
                      max=row$end,
                      value=c(row$start, row$end))
  }, ignoreInit = T)

  observeEvent(input$update, {
    output$geneModelPlot <- renderPlot({
      make_figure(ichrom=isolate(input$chrom), 
                  igene=isolate(input$gene), 
                  left=isolate(input$position[1]), 
                  right=isolate(input$position[2])
                  )
      })
    })
}




# RUN ---------------------------------------------------------------------

shinyApp(ui, server)


