library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(plotly)


targets <- read_tsv("targets.tsv")

ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq gene expression", titleWidth = 275),
  dashboardSidebar(
    width = 275,
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Quality control", tabName = "qc"),
      menuItem("Normalization", tabName = "norm"),
      menuItem("Results", tabName = "results")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "input",
        h1("Input"),
        p("TAC-seq data analysis output file with:"),
        tags$ul(
          tags$li("column sample"),
          tags$li("column locus"),
          tags$li("column molecule_count")
        ),
        fileInput("file", label = "Choose a file", accept = "text"),
        h3("Table"),
        dataTableOutput("table")
      ),
      tabItem(
        tabName = "qc",
        h1("Quality control"),
        p("ERCC spike-in controls"),
        h3("Plot"),
        plotOutput("spike_ins")
      ),
      tabItem(
        tabName = "norm",
        h1("Normalization"),
        p("Normalization by geometric mean of housekeeping genes"),
        h3("Plot"),
        plotOutput("housekeepers"),
        h3("Table"),
        dataTableOutput("norm_counts"),
        downloadButton("download", "Download normalized molecule counts")
      ),
      tabItem(
        tabName = "results",
        h1("Results"),
        p("Normalized molecule counts of targeted biomarkers"),
        h3("Plot"),
        plotlyOutput("pca")
      )
    )
  )
)

server <- function(input, output) {

  data <- reactive({
    req(input$file)
    read_tsv(input$file$datapath) %>%
      mutate(sample = str_remove(sample, "_S\\d+_*_R\\d_\\d+")) %>%
      filter(
        sample != "Undetermined",  # remove undetermined samples
        locus != "unmatched"  # remove unmatched loci
      ) %>%
      left_join(targets, by = c("locus" = "id"))
  })

  output$table <- renderDataTable({
    data()
  })

  output$spike_ins <- renderPlot({
    data() %>%
      filter(type == "spike_in") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      labs(title = "ERCC spike-ins", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  })

  output$housekeepers <- renderPlot({
    data() %>%
      filter(type == "housekeeper") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      labs(title = "Housekeeping genes", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  })

  norm_counts <- reactive({
    # normalize molecule counts -----------------------------------------------

    lst <- data() %>%
      select(sample, locus, molecule_count, type) %>%
      spread(sample, molecule_count) %>%
      split(.$type) %>%
      map(select, -type) %>%
      map(column_to_rownames, "locus") %>%
      map(as.matrix)

    bm <- lst$biomarker
    hk <- lst$housekeeper

    hk_geo_means <- exp(apply(log(hk), 2, mean))
    norm_counts <- sweep(bm, 2, hk_geo_means, "/")

    # remove housekeeper outliers ---------------------------------------------

    hk_zeros <- which(hk_geo_means == 0)
    if (length(hk_zeros) > 0) {
      norm_counts <- norm_counts[, -hk_zeros]
      cat("removed sample(s) where geometric mean of housekeeping genes is zero:", names(hk_zeros), sep = "\n")
    }

    t(norm_counts)
  })

  output$norm_counts <- renderDataTable({
    norm_counts()
  })

  output$pca <- renderPlotly({
    pca_norm <- norm_counts() %>%
      prcomp(scale. = T)

    pca_norm$x %>%
      as_tibble(rownames = "sample") %>%
      ggplot(aes(PC1, PC2, color = sample)) +
      geom_point() +
      labs(title = "Biomarkers", subtitle = "normalized molecule counts")
  })

  output$download <- downloadHandler(
    filename = "norm_counts.tsv",
    content = function(file) {
      norm_counts() %>%
      as_tibble(rownames = "sample") %>%
      write_tsv(file)
    }
  )

}

shinyApp(ui = ui, server = server)
