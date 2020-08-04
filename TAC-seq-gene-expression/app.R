library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(recipes)
library(plotly)
library(embed)
library(pheatmap)


ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq gene expression"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Quality control", tabName = "qc"),
      menuItem("Normalization", tabName = "norm"),
      menuItem("Visualization", tabName = "visual"),
      tags$a(
        href = "https://github.com/seqinfo/tac-seq-gene-expression", "GitHub",
        align = "center", style = "
        position:absolute;
        bottom:0;
        width:100%;
        weight:40px;   /* Height of the footer */
        color: white;
        padding: 10px;
        background-color: black;
        z-index: 1000;"
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "input",
        h1("Input"),
        p("Use TAC-seq data analysis output file as an input."),
        fileInput("counts", label = "Choose input file(s):", multiple = TRUE,
          accept = "text"
        ),
        selectInput(
          "target_list", label = "Choose target list or file:",
          choices = list(
            "Choose one" = "",
            "READY65 targets" = "data/targets/READY65_targets.tsv",
            "READY76 targets" = "data/targets/READY76_targets.tsv"
          )
        ),
        fileInput("target_file", label = "", accept = "text"),
        selectInput(
          "control_list", label = "Choose control list or file (optional):",
          choices = list(
            "Choose one" = "",
            "READY65 control set" = "data/controls/READY65_control_set.tsv",
            "READY65 small control set" = "data/controls/READY65_small_control_set.tsv"
          )
        ),
        fileInput("control_file", label = "", accept = "text"),
        plotOutput("biomarkers"),
        h2("Counts"),
        dataTableOutput("counts"),
        h2("Targets"),
        tableOutput("targets"),
        h2("Controls"),
        dataTableOutput("controls")
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
        h2("Normalized counts"),
        p("Samples, which geometric mean of housekeeping genes is zero, are removed."),
        dataTableOutput("norm_counts"),
        uiOutput("download")
      ),
      tabItem(
        tabName = "visual",
        h1("Visualization"),
        p("Visualizing the normalized molecule counts of targeted biomarkers."),
        plotlyOutput("pca", width = "auto", height = 700),
        plotlyOutput("umap", width = "auto", height = 700),
        plotOutput("heatmap", width = "auto", height = 700)
      )
    )
  )
)

server <- function(input, output) {

  geo_mean <- function(x) {
    # Compute the sample geometric mean.
    exp(mean(log(x)))
  }


# counts ------------------------------------------------------------------

  counts <- reactive({
    req(input$counts)

    validate(
      need(
        try(
          counts <- input$counts$datapath %>%
            set_names(nm = input$counts$name) %>%
            map_dfr(read_tsv, .id = "file") %>%
            select(file, sample, locus, molecule_count) %>%
            filter(
              !str_detect(sample, "Undetermined"),  # remove undetermined samples
              locus != "unmatched"  # remove unmatched loci
            )
        ),
        "Incorrect input file(s). Please choose correct TAC-seq count file(s) with columns 'sample', 'locus', and 'molecule_count'."
      )
    )

    validate(
      need(
        counts %>%
          count(sample, locus) %>%
          pull(n) == 1,
        str_c(
          "Remove duplicated sample: ",
          counts %>%
            count(sample, locus) %>%
            filter(n > 1) %>%
            pull(sample) %>%
            unique()
        )
      )
    )

    counts
  })

  output$counts <- renderDataTable(counts())


# targets -----------------------------------------------------------------

  targets <- reactive({
    if (isTruthy(input$target_file)) {
      validate(
        need(
          try(
            targets <- read_tsv(input$target_file$datapath) %>%
              select(target, type)
          ),
          "Incorrect target file. Please choose correct target file with columns 'target', and 'type'."
        )
      )
    } else if (isTruthy(input$target_list)) {
      targets <- read_tsv(input$target_list)
    } else {
      req(FALSE)
    }

    targets
  })

  output$targets <- renderTable(
    targets() %>%
      count(type)
  )


# controls ----------------------------------------------------------------

  controls <- reactive({
    biomarkers <- targets() %>%
      filter(type == "biomarker") %>%
      pull(target)

    if (isTruthy(input$control_file)) {
      validate(
        need(
          try(
            controls <- read_tsv(input$control_file$datapath) %>%
              select(sample, label, !!biomarkers)
          ),
          "Incorrect control file. Please choose correct control file with columns 'sample', 'label', and column for each 'target'."
        )
      )
    } else if (isTruthy(input$control_list)) {
      validate(
        need(
          try(
            controls <- read_tsv(input$control_list) %>%
              select(sample, label, !!biomarkers)
          ),
          "Missing target(s) in controls. Please choose controls with correct targets."
        )
      )
    } else {
      return(NULL)
    }

    controls
  })

  output$controls <- renderDataTable(controls())


# raw counts --------------------------------------------------------------

  raw_counts <- reactive(
    counts() %>%
      left_join(targets(), by = c("locus" = "target"))
  )

  output$raw_counts <- renderDataTable(raw_counts())


# biomarkers --------------------------------------------------------------

  output$biomarkers <- renderPlot(
    raw_counts() %>%
      filter(type == "biomarker") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      scale_y_log10() +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "Biomarkers", x = NULL, y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# spike-ins ---------------------------------------------------------------

  output$spike_ins <- renderPlot(
    raw_counts() %>%
      filter(type == "spike_in") %>%
      ggplot(aes(sample, molecule_count, color = locus)) +
      geom_point() +
      geom_line(aes(group = locus)) +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "ERCC spike-in controls", x = NULL, y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# housekeepers ------------------------------------------------------------

  output$housekeepers <- renderPlot(
    raw_counts() %>%
      filter(type == "housekeeper") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_point(aes(color = locus)) +
      stat_summary(aes(statistic = "geometric mean"), fun = geo_mean,
                   geom = "crossbar") +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "Housekeeping genes", x = NULL, y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# normalization -----------------------------------------------------------

  norm_counts <- reactive(
    raw_counts() %>%
      group_by(sample) %>%
      mutate(norm_factor = geo_mean(molecule_count[type == "housekeeper"])) %>%
      ungroup() %>%
      filter(type == "biomarker") %>%
      mutate(norm_molecule_count = molecule_count / norm_factor) %>%
      filter(is.finite(norm_molecule_count)) %>%
      pivot_wider(id_cols = sample, names_from = locus,
                  values_from = norm_molecule_count)
  )

  output$norm_counts <- renderDataTable(norm_counts())


# download ----------------------------------------------------------------

  output$download <- renderUI({
    req(norm_counts())
    downloadButton("file", "Download normalized counts")
  })

  output$file <- downloadHandler(
    filename = "TAC-seq_normalized_counts.tsv",
    content = function(file) {
      norm_counts() %>%
        write_tsv(file)
    }
  )


# train and test data -----------------------------------------------------

  train_data <- reactive(req(controls()))
  test_data <- reactive(req(bind_rows(controls(), norm_counts())))


# PCA ---------------------------------------------------------------------

  output$pca <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 2) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(PC1, PC2, color = label, sample = sample)) +
      geom_point() +
      labs(title = "PCA of biomarkers", color = NULL) +
      coord_equal()
    ggplotly()
  })


# UMAP --------------------------------------------------------------------

  output$umap <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_umap(all_numeric(), seed = c(1, 1)) %>%
      # step_string2factor(label) %>%
      # step_umap(all_numeric(), outcome = vars(label), seed = c(1, 1)) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(umap_1, umap_2, color = label, sample = sample)) +
      geom_point() +
      labs(title = "UMAP of biomarkers") +
      coord_equal()
    ggplotly()
  })


# heatmap -----------------------------------------------------------------

  output$heatmap <- renderPlot(
    test_data() %>%
      column_to_rownames("sample") %>%
      select(where(is.numeric)) %>%
      t() %>%
      na_if(0) %>%
      log() %>%
      pheatmap(treeheight_row = 0, main = "Biomarkers")
  )
}


shinyApp(ui = ui, server = server)
