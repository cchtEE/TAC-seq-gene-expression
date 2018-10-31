library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(plotly)
library(pheatmap)


targets <- read_tsv("data/targets/READY_targets.tsv")
reference <- read_tsv("data/controls/READY_controls.tsv")

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
        p("Use TAC-seq data analysis output file as an input."),
        fileInput("file", label = "Choose a file", accept = "text"),
        plotOutput("biomarkers"),
        h2("Raw data"),
        dataTableOutput("table")
      ),
      tabItem(
        tabName = "qc",
        h1("Quality control"),
        p("ERCC spike-in controls are used for quality control."),
        plotOutput("spike_ins")
      ),
      tabItem(
        tabName = "norm",
        h1("Normalization"),
        p("Molecule counts are normalized by geometric mean of housekeeping genes."),
        plotOutput("housekeepers"),
        h2("Normalized data"),
        p("Samples, which geometric mean of housekeeping genes is zero, are removed."),
        dataTableOutput("bm_norm"),
        downloadButton("download", "Download normalized data")
      ),
      tabItem(
        tabName = "results",
        h1("Results"),
        p("Visualizing the normalized molecule counts of targeted biomarkers."),
        plotlyOutput("pca")
      )
    )
  )
)

server <- function(input, output) {

  data <- reactive({
    req(input$file)
    read_tsv(input$file$datapath) %>%
      filter(
        !str_detect(sample, "Undetermined"),  # remove undetermined samples
        locus != "unmatched"  # remove unmatched loci
      ) %>%
      left_join(targets, by = c("locus" = "id"))
  })

  output$table <- renderDataTable({
    data()
  })

  output$biomarkers <- renderPlot({
    data() %>%
      filter(type == "biomarker") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      scale_y_log10() +
      labs(title = "Biomarkers", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  })

  output$spike_ins <- renderPlot({
    data() %>%
      filter(type == "spike_in") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      labs(title = "ERCC spike-in controls", subtitle = "raw molecule counts", y = "molecule count") +
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

  bm_norm <- reactive({
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
    bm_norm <- sweep(bm, 2, hk_geo_means, "/")

    # remove housekeeper outliers ---------------------------------------------

    hk_zeros <- which(hk_geo_means == 0)
    if (length(hk_zeros) > 0) {
      bm_norm <- bm_norm[, -hk_zeros]
    }

    t(bm_norm)
  })

  output$bm_norm <- renderDataTable({
    bm_norm()
  })

  output$download <- downloadHandler(
    filename = "normalized_counts.tsv",
    content = function(file) {
      bm_norm() %>%
        as_tibble(rownames = "sample") %>%
        write_tsv(file)
    }
  )

  output$pca <- renderPlotly({
    pca_norm <- bm_norm() %>%
      prcomp(scale. = T)

    pca_norm$x %>%
      as_tibble(rownames = "sample") %>%
      right_join()
      ggplot(aes(PC1, PC2, color = sample)) +
      geom_point() +
      labs(title = "Biomarkers")
  })

}

shinyApp(ui = ui, server = server)
