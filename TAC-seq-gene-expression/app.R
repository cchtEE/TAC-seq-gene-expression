library(shiny)
library(shinydashboard)
library(tidyverse)
library(DT)
library(recipes)
library(pheatmap)
library(plotly)
library(embed)


ui <- dashboardPage(
  dashboardHeader(title = "TAC-seq gene expression"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input"),
      menuItem("Normalization", tabName = "normalization"),
      menuItem("Visualization", tabName = "visualization"),
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
            "READY61 targets" = "data/targets/READY61_targets.tsv",
            "READY72 targets" = "data/targets/READY72_targets.tsv"
          )
        ),
        fileInput("target_file", label = "", accept = "text"),
        h3("Counts"),
        plotOutput("count_plot"),
        dataTableOutput("counts"),
        h3("Targets"),
        dataTableOutput("targets")
      ),
      tabItem(
        tabName = "normalization",
        h1("Normalization"),
        p("Molecule counts are normalized by geometric mean of housekeeping genes."),
        h3("Housekeeping genes"),
        plotOutput("housekeepers"),
        h3("Normalized counts"),
        p("Samples with geometric mean of zero are removed."),
        dataTableOutput("norm_biomarkers"),
        uiOutput("download")
      ),
      tabItem(
        tabName = "visualization",
        h1("Visualization"),
        p("Visualizing the normalized molecule counts of targeted biomarkers."),
        selectInput(
          "control_list", label = "Choose control list or file (optional):",
          choices = list(
            "Choose one" = "",
            "READY65 control set" = "data/controls/READY65_control_set.tsv",
            "READY65 small control set" = "data/controls/READY65_small_control_set.tsv"
          )
        ),
        fileInput("control_file", label = "", accept = "text"),
        h3("Controls"),
        dataTableOutput("controls"),
        h3("Heatmap"),
        plotOutput("heatmap", width = "auto", height = 700),
        h3("PCA"),
        plotlyOutput("pca", width = "auto", height = 700),
        h3("UMAP"),
        plotlyOutput("umap", width = "auto", height = 700)
      )
    )
  )
)

server <- function(input, output) {


# counts ------------------------------------------------------------------

  counts <- reactive({
    count_file <- input$counts

    req(count_file, targets())

    n_targets = nrow(targets())

    validate(
      need(
        try(
          counts <- count_file$datapath %>%
            set_names(nm = count_file$name) %>%
            map_dfr(read_tsv, .id = "file") %>%
            filter(!str_detect(sample, "Undetermined"),  # remove undetermined samples
                   locus != "unmatched") %>%  # remove unmatched loci
            right_join(targets(), by = c("locus" = "target")) %>%
            complete(nesting(file, sample), nesting(locus, type)) %>%
            group_by(sample) %>%
            # filter(n() == n_targets,  # remove duplicated samples
            #        !any(is.na(molecule_count))) %>%  # remove samples with missing targets
            mutate(hk_geo_mean = exp(mean(log(molecule_count[type == "housekeeper"]))),  # geometric mean of housekeeping genes
                   norm_molecule_count = molecule_count / hk_geo_mean) %>%
            ungroup()
        ),
        "Incorrect input file(s). Please choose correct TAC-seq count file(s) with columns 'sample', 'locus', and 'molecule_count'."
      )
    )

    validate(
      need(
        counts %>%
          count(sample) %>%
          filter(n != n_targets) %>%
          nrow() == 0,
        counts %>%
          count(sample) %>%
          filter(n != n_targets) %>%
          transmute(duplicated = str_c(sample, " is duplicated")) %>%
          pull()
      ),
      need(
        !anyNA(counts$molecule_count),
        counts %>%
          filter(is.na(molecule_count), !is.na(sample)) %>%
          transmute(missing = str_c(sample, ": ", locus, " is missing")) %>%
          pull()
      )
    )

    counts
  })

  output$counts <- renderDataTable(counts())


# targets -----------------------------------------------------------------

  targets <- reactive({
    target_file <- input$target_file
    target_list <- input$target_list

    if (isTruthy(target_file)) {
      validate(
        need(
          try(
            targets <- read_tsv(target_file$datapath) %>%
              select(target, type)
          ),
          "Incorrect target file. Please choose correct target file with columns 'target', and 'type'."
        )
      )
    } else if (isTruthy(target_list)) {
      targets <- read_tsv(target_list)
    } else {
      req(FALSE)
    }

    targets
  })

  output$targets <- renderDataTable(
    targets()
  )


# plot counts -------------------------------------------------------------

  output$count_plot <- renderPlot(
    counts() %>%
      ggplot(aes(sample, molecule_count)) +
      geom_boxplot() +
      geom_point(aes(color = locus), show.legend = FALSE) +
      scale_y_log10() +
      facet_wrap(vars(file), scales = "free") +
      labs(x = NULL, y = "molecule count", color = NULL) +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# plot housekeepers -------------------------------------------------------

  output$housekeepers <- renderPlot(
    counts() %>%
      filter(type == "housekeeper") %>%
      ggplot(aes(sample, molecule_count)) +
      geom_point(aes(color = locus)) +
      geom_errorbar(aes(y = hk_geo_mean, ymin = hk_geo_mean, ymax = hk_geo_mean),
                    size = 1) +
      facet_wrap(vars(file), scales = "free") +
      labs(x = NULL, y = "molecule count") +
      theme(axis.text.x = element_text(vjust = 0.5, angle = 90))
  )


# biomarker counts --------------------------------------------------------

  norm_biomarkers <- reactive({
    counts() %>%
      filter(type == "biomarker") %>%
      pivot_wider(id_cols = c(file, sample), names_from = locus,
                  values_from = norm_molecule_count) %>%
      filter(across(where(is.numeric), is.finite))
  })

  output$norm_biomarkers <- renderDataTable(
    norm_biomarkers(),
    options = list(scrollX = TRUE)
  )


# download ----------------------------------------------------------------

  output$download <- renderUI({
    req(norm_biomarkers())
    downloadButton("file", "Download normalized counts")
  })

  output$file <- downloadHandler(
    filename = "TAC-seq_normalized_counts.tsv",
    content = function(file) {
      norm_biomarkers() %>%
        write_tsv(file)
    }
  )


# controls ----------------------------------------------------------------

  controls <- reactive({
    control_file <- input$control_file
    control_list <- input$control_list

    req(targets())

    biomarkers <- targets() %>%
      filter(type == "biomarker") %>%
      pull(target)

    if (isTruthy(control_file)) {
      validate(
        need(
          try(
            controls <- read_tsv(control_file$datapath) %>%
              select(sample, group, !!biomarkers)
          ),
          "Incorrect control file. Please choose correct control file with columns 'sample', 'group', and column for each 'target'."
        )
      )
    } else if (isTruthy(control_list)) {
      validate(
        need(
          try(
            controls <- read_tsv(control_list) %>%
              select(sample, group, !!biomarkers)
          ),
          "Missing target(s) in controls. Please choose correct controls or targets."
        )
      )
    } else {
      return(NULL)
    }

    controls
  })

  output$controls <- renderDataTable(controls(), options = list(scrollX = TRUE))


# train and test data -----------------------------------------------------

  train_data <- reactive(req(controls()))
  test_data <- reactive({
    req(controls())

    bind_rows(norm_biomarkers(), controls())
  })


# heatmap -----------------------------------------------------------------

  output$heatmap <- renderPlot(
    test_data() %>%
      column_to_rownames("sample") %>%
      select(where(is.numeric)) %>%
      t() %>%
      na_if(0) %>%
      log() %>%
      pheatmap(treeheight_row = 0)
  )


# PCA ---------------------------------------------------------------------

  output$pca <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      step_pca(all_numeric(), num_comp = 2) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(PC1, PC2, color = group, sample = sample)) +
      geom_point() +
      coord_equal()
    ggplotly()
  })


# UMAP --------------------------------------------------------------------

  output$umap <- renderPlotly({
    train_data() %>%
      recipe() %>%
      step_normalize(all_numeric()) %>%
      # step_umap(all_numeric(), seed = c(1, 1)) %>%
      step_string2factor(group) %>%
      step_umap(all_numeric(), outcome = vars(group), seed = c(1, 1)) %>%
      prep(strings_as_factors = FALSE) %>%
      bake(new_data = test_data()) %>%
      ggplot(aes(umap_1, umap_2, color = group, sample = sample)) +
      geom_point() +
      coord_equal()
    ggplotly()
  })
}


shinyApp(ui = ui, server = server)
