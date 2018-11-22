library(shiny)
library(shinydashboard)
library(DT)
library(tidyverse)
library(plotly)
library(pheatmap)


ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq gene expression", titleWidth = 275),
  dashboardSidebar(
    width = 275,
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Quality control", tabName = "qc"),
      menuItem("Normalization", tabName = "norm"),
      menuItem("Visualization", tabName = "visual"),
      tags$a(href = "https://github.com/cchtEE/TAC-seq-gene-expression", "GitHub", align = "center", style = "
              position:absolute;
              bottom:0;
              width:100%;
              height:40px;   /* Height of the footer */
              color: white;
              padding: 10px;
              background-color: black;
              z-index: 1000;")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "input",
        h1("Input"),
        p("Use TAC-seq data analysis output file as an input."),
        fileInput("counts", label = "Choose count file(s):", accept = "text", multiple = T),
        selectInput("target_list", label = "Choose target list or file:", choices = list(
          "Choose one" = "",
          "READY targets" = "data/targets/READY_targets.tsv"
        )),
        fileInput("target_file", label = "", accept = "text"),
        selectInput("control_list", label = "Choose control list or file (optional):", choices = list(
          "Choose one" = "",
          "READY controls" = "data/controls/READY_controls.tsv"
        )),
        fileInput("control_file", label = "", accept = "text"),
        plotOutput("biomarkers"),
        h2("Targets"),
        tableOutput("targets"),
        h2("Raw data"),
        dataTableOutput("data"),
        h2("Control data"),
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
        h2("Normalized data"),
        p("Samples, which geometric mean of housekeeping genes is zero, are removed."),
        p("Removed sample(s):"),
        verbatimTextOutput("outliers"),
        dataTableOutput("bm"),
        uiOutput("download")
      ),
      tabItem(
        tabName = "visual",
        h1("Visualization"),
        p("Visualizing the normalized molecule counts of targeted biomarkers."),
        plotlyOutput("pca"),
        plotOutput("heatmap")
      )
    )
  )
)

server <- function(input, output) {

# counts ------------------------------------------------------------------

  counts <- reactive({
    req(input$counts)
    validate(need(try(
      counts <- input$counts$datapath %>%
        set_names(nm = input$counts$name) %>%
        map_dfr(read_tsv, .id = "file") %>%
        select(file, sample, locus, molecule_count) %>%
        filter(
          !str_detect(sample, "Undetermined"),  # remove undetermined samples
          locus != "unmatched"  # remove unmatched loci
        )
      ),
      "Incorrect count file(s). Please choose correct TAC-seq count file(s) with columns \"sample\", \"locus\" and \"molecule_count\"."
    ))
    counts
  })

  output$counts <- renderDataTable(counts())

# targets -----------------------------------------------------------------

  targets <- reactive({
    if (isTruthy(input$target_file)) {
      validate(need(
        try(
          targets <- input$target_file$datapath %>%
            read_tsv() %>%
            select(target, type)
        ),
        "Incorrect target file. Please choose correct target file with columns \"target\" and \"type\"."
      ))
      targets
    } else if (isTruthy(input$target_list)) {
      input$target_list %>%
        read_tsv()
    } else {
      req(FALSE)
    }
  })

  output$targets <- renderTable(
    targets() %>%
      count(type)
  )

# controls ----------------------------------------------------------------

  controls <- reactive({
    bm <- targets() %>%
      filter(type == "biomarker") %>%
      pull(target)

    if (isTruthy(input$control_file)) {
      validate(need(
        try(
          controls <- input$control_file$datapath %>%
            read_tsv() %>%
            select(sample, label, !!bm)
        ),
        "Incorrect control file. Please choose correct control file with columns \"sample\", \"label\" and column for each \"target\"."
      ))
      controls
    } else if (isTruthy(input$control_list)) {
      input$control_list %>%
        read_tsv() %>%
        select(sample, label, !!bm)
    } else {
      return(NULL)
    }
  })

  output$controls <- renderDataTable(controls())

# data --------------------------------------------------------------------

  data <- reactive(
    counts() %>%
      left_join(targets(), by = c("locus" = "target"))
  )

  output$data <- renderDataTable(data())

# biomarkers --------------------------------------------------------------

  output$biomarkers <- renderPlot(
    data() %>%
      filter(type == "biomarker") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      scale_y_log10() +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "Biomarkers", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )

# spike-ins ---------------------------------------------------------------

  output$spike_ins <- renderPlot(
    data() %>%
      filter(type == "spike_in") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "ERCC spike-in controls", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )

# housekeepers ------------------------------------------------------------

  output$housekeepers <- renderPlot(
    data() %>%
      filter(type == "housekeeper") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus)) +
      facet_wrap(vars(file), scales = "free") +
      labs(title = "Housekeeping genes", subtitle = "raw molecule counts", y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )

# normalization -----------------------------------------------------------

  bm <- reactive({
    if (isTruthy(data()))
    validate(need(
      try(
        lst <- data() %>%
          select(sample, locus, molecule_count, type) %>%
          spread(sample, molecule_count) %>%
          split(.$type) %>%
          map(select, -type) %>%
          map(column_to_rownames, "locus") %>%
          map(as.matrix)
      ),
      "Duplicated samples. Please rename or remove duplicated samples from count files."
    ))

    bm_raw <- lst$biomarker
    hk_raw <- lst$housekeeper

    norm_factor <- exp(apply(log(hk_raw), 2, mean))  # housekeepers geometric mean
    bm <- sweep(bm_raw, 2, norm_factor, "/")

    # remove housekeeper outliers
    outliers <- "none"
    hk_zeros <- which(norm_factor == 0)
    if (length(hk_zeros) > 0) {
      outliers <- names(hk_zeros)
      bm <- bm[, -hk_zeros]
    }

    output$outliers <- renderPrint({
        cat(outliers, sep = "\n")
      })

    t(bm)
  })

  output$bm <- renderDataTable({
    bm()
  })

# download ----------------------------------------------------------------

  output$download <- renderUI({
    req(bm())
    downloadButton("file", "Download normalized data")
  })

  output$file <- downloadHandler(
    filename = "TAC-seq_normalized_counts.tsv",
    content = function(file) {
      bm() %>%
        as_tibble(rownames = "sample") %>%
        write_tsv(file)
    }
  )

# PCA ---------------------------------------------------------------------

  output$pca <- renderPlotly({
    if (isTruthy(controls())) {
      df <- bm() %>%
        as_tibble(rownames = "sample") %>%
        bind_rows(controls())

      mat <- df %>%
        select(-label) %>%
        column_to_rownames("sample") %>%
        as.matrix()

      pca <- mat %>%
        prcomp(scale. = T)

      pca$x %>%
        as_tibble(rownames = "sample") %>%
        left_join(df) %>%
        ggplot(aes(PC1, PC2, color = label, group = sample)) +
        geom_point() +
        labs(title = "Biomarkers")
    } else {
      pca <- bm() %>%
        prcomp(scale. = T)

      pca$x %>%
        as_tibble(rownames = "sample") %>%
        ggplot(aes(PC1, PC2, color = sample)) +
        geom_point() +
        labs(title = "Biomarkers")
    }
    ggplotly()
  })

# heatmap -----------------------------------------------------------------

  output$heatmap <- renderPlot({
    if (isTruthy(controls())) {
      mat <- bm() %>%
        as_tibble(rownames = "sample") %>%
        bind_rows(controls()) %>%
        select(-label) %>%
        column_to_rownames("sample") %>%
        as.matrix()
    } else {
      mat <- bm()
    }
    mat[mat == 0] <- NA  # replace 0 with NA
    pheatmap(log(t(mat)), scale = "row", treeheight_row = 0)
  })

}

shinyApp(ui = ui, server = server)
