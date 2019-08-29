library(shiny)
library(shinydashboard)
library(NMF)
library(shinyBS)
library(shinythemes)
nmf.options(grid.patch=TRUE)
#source("/Users/eqiao14/Desktop/Visualization/helpers.R")

# beginning of fluid page
dashboardPage(skin = "green",
  dashboardHeader(title = "Heatmap"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Getting Started",tabName = "tab2", icon = icon("coffee")
      ),
      menuItem("Heatmap",tabName = "tab1", icon = icon("sliders")),
      menuItem(
        "Upload Files", icon = icon('upload'),
        
        fileInput('file1', 'Upload Gene Table',
                  accept = c(
                    '.csv'
                  )
        ),
        fileInput('annotate_1', 'Upload Column Annotation', accept = c(".csv")),
        fileInput('annotate_4', 'Upload Row Annotation', accept = c(".csv"))),
      
    menuItem("Thresholds", icon = icon('scissors'),
             checkboxInput("Lower", label = "Lower Threshold", value = FALSE),
             checkboxInput("Upper", label = "Upper Threshold", value = FALSE),
             numericInput("Lower_T", label = h3("Lower Threshold"), value = -10),
             numericInput("Upper_T", label = h3("Upper Threshold"), value = 10)),
    menuItem("Normalization", icon = icon('magic'),
             checkboxInput("Normalize", label = "Normalize to HK", value = FALSE),
             checkboxInput("Normalize_med", label = "Normalize to Median Values", value = FALSE)),
    menuItem("Download", icon = icon('download'),
             textInput("text", label = ("Download Heatmap File Name"), value = "Enter name..."),
             downloadButton('downloadData', 'Download'),
             textInput("csv", label = ("Download Data Table File Name"), value = "Enter name..."),
             downloadButton('downloadData_2', 'Download')),
    menuItem("About",tabName = "tab3", icon = icon("info-circle")
    )
  )
),
    dashboardBody(
      tabItems(
        tabItem(tabName = "tab1",
          fluidPage(
            column(width = 7, class = "well",
             fluidRow(
               h4("Heatmap"),
               plotOutput("plot_b",  brush = brushOpts(
                 id = "plot2_brush"), height = 500, width = 500)
               ), bsAlert("no_HK")),
          column(width = 12, class = "well",
               fluidRow(
               h4("Zoomed Plot"),
               uiOutput('zoom')), 
               fluidRow(
               div(style = "height: 27px;",
                      sliderInput("Zoom", min = 500,
                                  max = 950, value = 950, label = NULL)))
          )
          )
        ),
        tabItem(tabName = 'tab2',
                
                h2("Data Visualization Tutorial", align = 'center'),
                h4("Step 1. Uploading Files"),
                fluidRow(
                  box(
                    title = "Gene table",
                    img(src = "Genes.png", height = 400, width = 400),
                    br(),
                    br(),
                    p("The Gene Table should have the gene names in first column followed by all of the samples.
                      These should all be .csv files.",
                      style = ""),
                    br()),
                    box(
                    title = "Annotation Files",
                    img(src = "Annotations.png", height = 400, width = 200),
                    br(),
                    br(),
                    p("The visualizer is able to generate annotations of each sample or transcript. These are uploaded
                      via the Column Annotation or Row Annotation options. Both follow the same format, with annotation
                       categories listed as columns within a .csv file.", style = ""))
                    ),
                br(),
                h4("Step 2. Heatmap Manipulation"),
                fluidRow(
                  box(
                    title = "Heatmap Output",
                    p("After uploading files, a heatmap should be generated similar to this one.", style = ""),
                    br(),
                    br(),
                    img(src = "Heatmap.png", height = 400, width = 400),
                    br(),
                    br(),
                    p('The column annotations show the sex of each sample, the row annotations show the gene function.', 
                      style = "")
                  ),
                  box(
                    title = "Heatmap Thresholds",
                    p("Upper and lower thresholds can be set, which set values to respective thresholds. This shows
                      the same data with the original heatmap but with an lower threshold expression of -10 and upper threshold
                      expression of 10.", style = ''),
                    img(src = "Heatmap_thresholds.png", height = 400, width = 400),
                    br(),
                    br(),
                    p('To select thresholds, select the Thresholds option on the sidebar and input thresholds. 
                      Check the corresponding threshold checkbox to see changes to heatmap.',
                      style = '')
                    
                  ),
                  box(
                    title = "Median Centering",
                    p("Two options exist to normalize data, median centering and normalizing to a
                      specific gene included in the expression data. Median centering will subtract
                      the median value from each transcript. Using the same data, the following heatmap
                      was generated.", style = ''),
                    img(src = "Heatmap_median.png", height = 400, width = 400),
                    br(),
                    br(),
                    p('Select the \'Normalize to Median Values\' option if median centering is required.',
                      style = ''),
                    br(),
                    br()
                    ),
                  box(
                    title = "Normalizing to Gene",
                    p('To normalize to a gene in the data, transcript annotation is required. Label any
                      gene \'HK\' that is intended to be used as a housekeeping gene or gene for normalizing.
                      These labels should be in the first column of the .csv file for row annotations.
                      ', style = ''),
                    img(src = "Heatmap_HK.png", height = 400, width = 400),
                    br(),
                    br(),
                    p('Select the \'Normalize to HK\' option if such normalization is required.', style = '
                      '),
                    p('Note: Normalization
                      is based on log2 of expression values. Using raw read count 
                      will not have the intended effect.', style = '')
                    )
                ),
                h4('Step 3. Downloading Heatmaps'),
                fluidRow(
                  box(
                    title = 'Formats', width = 800,
                    p('Heatmaps can be downloaded as a PDF. Click the Download option, name the files and click the download button. 
                      In addition to these, the values used to generate the 
                      heatmap can also be downloaded as a csv file. For instance, if a threshold was used, 
                      a table with the manipulated values can be downloaded.', style = ''),
                    p('The sample data used to generate all these plots can be found', a('here', href=
                                                                                         'https://www.mediafire.com/?7d5uklmf48ml5oo',
                                                                                         target="_blank"))
                      
                  )
                )
                
        ),
        tabItem(tabName = 'tab3',
          box(title = 'References', width = 800,
              p('The algorithm for clustering is based on the NMF package in R.'),
              h6('Renaud Gaujoux, Cathal Seoighe (2010). A flexible R package for nonnegative matrix factorization. BMC
                 Bioinformatics 2010, 11:367. [http://www.biomedcentral.com/1471-2105/11/367]'))
        )
        )
      )
    )





