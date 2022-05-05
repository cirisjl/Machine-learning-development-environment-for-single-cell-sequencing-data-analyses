# ----------------------------------------
# --       PROGRAM ui_sidebar.R         --
# ----------------------------------------
# USE: Create UI elements for the
#      application sidebar (left side on
#      the desktop; contains options) and
#      ATTACH them to the UI by calling
#      add_ui_sidebar_basic() or
#      add_ui_sidebar_advanced()
#
# NOTEs:
#   - All variables/functions here are
#     not available to the UI or Server
#     scopes - this is isolated
# ----------------------------------------

# -- IMPORTS --
suppressPackageStartupMessages(library(shinyjs))


# ----------------------------------------
# --     SIDEBAR ELEMENT CREATION       --
# ----------------------------------------

# -- Create Basic Elements

sidebar_header <- tags$div(useShinyjs())

file_input <- tags$div(
    style = "margin-top:20px;margin-bottom:25px;",
    align = "center",
    fileInput(inputId = "fileInputDialog",
              label   = NULL,
              buttonLabel = "Load New Data Object",
              width   = "80%"))

data_select <- tags$div(
  style = "margin-bottom:20px;",
  selectInput(inputId = "dataset", "Choose dataset", choices=c("Allen H. MTG", "Allen H. VC", "Allen H. ACC", "Leng et al. H. EC_Exc", "Leng et al. H. EC_Inh", "Leng et al. H. SFG_Exc", "Leng et al. H. SFG_Inh",
                                                               "Allen M. VC", "Allen M. ALM", "Allen M. ACA", 
                                                               "Allen M. MOpC", "Allen M. MOpN")))
# selectInput(inputId = "dataset", "Choose dataset", choices=c("Allen Human Middle Temporal Gyrus", "Allen Human Primary Visual Cortex", "Allen Human Anterior Cingulate Cortex",
#"Allen Mouse Primary Visual Cortex", "Allen Mouse Anterior Lateral Motor Area", "Allen Mouse Anterior Cingulate Area", 
#"Allen Mouse Primary Motor Area cells", "Allen Mouse Primary Motor Area nuclei", "SFG"))) "H. SFG", "H. EC", 

genes_select <- tags$div(
    #h4("Chart Options"),
    selectizeInput(inputId = "gene",
                   label    = "Enter your target",
                   choices  = NULL,
                   multiple = FALSE,
                   options  = list(maxOptions = 5,
                                   placeholder = "Type/Click then Select",
                                   searchField = "value",
                                   plugins     = list('remove_button'))))

help_text <- tags$div(
    tags$br(),
    tags$h4(tags$a(href = 'canvasxpress.org', "CanvasXpress"), "Tips"),
    tags$p(style = "margin:10px;",
           tags$strong('Toolbar'), "- the main toolbar is available if you hover over the top title of the chart. Additional functionality can also be accessed by right-clicking anywhere on the chart",
           tags$br(), tags$br(),
           tags$strong('Zoom'), "- select an area or use the mouse to zoom by scrolling.  Reset the canvas by hitting Esc",
           tags$br(), tags$br(),
           tags$strong('Focus'), "- select a legend item to toggle a fade on that item",
           tags$br(), tags$br(),
           tags$strong('Select'), "- select points on a plot using shift-drag to select an area.  Deselect by clicking a blank area or hitting Esc",
           tags$br(), tags$br(),
           tags$strong('Download'), "- hover over the main title and select the camera icon from the toolbar"))

about_text <- tags$div(
                tags$br(),
                tags$h4(style = "margin:-10px;",
                        actionLink("about_link", "About This App")))

#                          help_text,
add_ui_sidebar_basic(list(sidebar_header,
                         # file_input,
                          data_select,
                          genes_select,
                          about_text),
                     tabname = "Application")

# -- Create Advanced Elements

#add_ui_sidebar_advanced(list(uiOutput("filterOptions"),
  #                           hr()),
   #                     tabname = "Filtering")

