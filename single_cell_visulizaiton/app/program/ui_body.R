# ----------------------------------------
# --          PROGRAM ui_body.R         --
# ----------------------------------------
# USE: Create UI elements for the
#      application body (right side on the
#      desktop; contains output) and
#      ATTACH them to the UI by calling
#      add_ui_body()
#
# NOTEs:
#   - All variables/functions here are
#     not available to the UI or Server
#     scopes - this is isolated
# ----------------------------------------

# -- IMPORTS --



# ----------------------------------------
# --      BODY ELEMENT CREATION         --
# ----------------------------------------

# -- Create Elements

# Add css class to radiobuttons div to position them on same row as the text & handle behaviour related to the file upload.
application_js <- '$(document).on("shiny:connected", function(e) {
                    var elements = document.getElementsByClassName("shiny-options-group");
                    for (var i = 0; i < elements.length; i++)
                    {
                        elements[i].classList.add("radio-inline");
                    }
                   });
                   $(document).on("change", "#fileInputDialog", function(e) {
                        Shiny.onInputChange("fileChosen", this.files[0].size);
                   });
                   $(document).on("focusout", "#loading_modal", function(e) {
                        Shiny.onInputChange("fileModalClosed", Math.random());
                   });
                   Shiny.addCustomMessageHandler("openTitleInfoBox",
                        function(message) {
                           $("#titleinfobox_trigger").trigger("click");
                        }
                    );
                    
var dimension = [0, 0];
$(document).on("shiny:connected", function(e) {
  dimension[0] = window.innerWidth;
  dimension[1] = window.innerHeight;
  Shiny.onInputChange("dimension", dimension);
});
$(window).resize(function(e) {
  dimension[0] = window.innerWidth;
  dimension[1] = window.innerHeight;
  Shiny.onInputChange("dimension", dimension);
});'

loading_modal <- bsModal("loading_modal",
                           title   = "SCRSDE Loading...",
                           trigger = NULL,
                           size    = "small",
                           tags$head(tags$style("#loading_modal .modal-title {color:#3380A9; font-weight:bold;}"),
                                     tags$style("#loading_modal .modal-footer {display:none;}")),
                           h4(align = 'center', 'The file you have selected is larger than 10 MB. It is now being uploaded and will load momentarily.'),
                           br(),
                           p(align = 'center', em('Click anywhere outside this dialog to continue')))

file_error_modal <- bsModal("file_error_modal",
                            title   = "SCRSDE File Upload Error",
                            trigger = NULL,
                            size    = "large",
                            tags$head(tags$style("#file_error_modal .modal-title {color:#3380A9; font-weight:bold; text-align:center;}"),
                                      tags$style("#file_error_modal .modal-footer {display:none;}"),
                                      tags$style("#missing_list {width:60%; margin-left:25%}"),
                                      tags$style("#missing_list li {text-align:left;}")),
                            h4(align = 'center', htmlOutput("loading_error_message")),
                            br(),
                            p(align = 'center', em('Click anywhere outside this dialog to continue')))

create_bs_button <- function(id, label, width = "100%", disabled = FALSE) {
    bsButton(inputId  = id,
             label    = label,
             type     = "primary",
             value    = FALSE,
             style    = "primary",
             size     = "default",
             width    = width,
             icon     = icon("arrow-right"),
             block    = FALSE,
             disabled = disabled)
}

body1 <- box(id          = "bodyElement1",
             title       = "Cluster/Gene Expression Plots",
             width       = 14,
             collapsible = FALSE,
             tags$head(tags$script(HTML(application_js))),
             tags$link(rel = "stylesheet", type = "text/css", href = "style.css"),
             tabBox(id     = "tabBodyElement2",
                    width  = 14,
                    selected = "UMAP/tSNE",
                    tabPanel("UMAP/tSNE",
                             fluidRow(
                                 #column(width = 3, selectInput(inputId = "dataset", "Choose dataset", choices=c("Allen Human Middle Temporal Gyrus", "Allen Human Primary Visual Cortex", "Allen Human Anterior Cingulate Cortex",
                                  #                                                                              "Allen Mouse Primary Visual Cortex", "Allen Mouse Anterior Lateral Motor Area", "Allen Mouse Anterior Cingulate Area", 
                                   #                                                                             "Allen Mouse Primary Motor Area cells", "Allen Mouse Primary Motor Area nuclei"))),
                                 column(width = 3, radioButtons(inputId = "plot", label = "Plot type", choices = c("UMAP", "tSNE"))),
                                 column(width = 3, radioButtons(inputId = "dim", label = "Dimension", choices = c("2D", "3D"),inline = T)),
                                 column(width = 4, radioButtons(inputId = "lay", label = "Plot circles layout", choices = c("Clusters", "Gene Expression"))),
                                 column(width = 2, tags$label(HTML("&zwnj;")),
                                                    create_bs_button(id    = "go",
                                                                     label = " Plot/Refresh"))),
                             tags$hr(),
                             fluidRow(
                               column(width = 12, plotlyOutput("tsne.plot", height = "auto"))
                             ),
                             tags$hr(),
                             fluidRow(
                               column(width = 6, tableOutput("tsne.click"))
                             ))
                    # tabPanel("Compare clusters",
                    #          fluidRow(
                    #            column(width = 4, selectInput(inputId = "dataseto", "Compare with", choices=c("Leng et al. H. EC_Exc", "Leng et al. H. EC_Inh",
                    #                                                                                          "Leng et al. H. SFG_Exc", "Leng et al. H. SFG_Inh", "Allen H. VC", "Allen H. ACC", 
                    #                                                                                          "Allen M. VC", "Allen M. ALM", "Allen M. ACA", 
                    #                                                                                          "Allen M. MOpC", "Allen M. MOpN"))),
                    #            column(width = 3, tags$label(HTML("&zwnj;")), tags$br(), create_bs_button(id    = "heatmapPlotBtn",
                    #                                                 label = " Plot/Refresh"))),
                    #          tags$hr(),
                    #          fluidRow(
                    #            column(width = 12, plotlyOutput(outputId = 'heatmap', height = "auto"))
                    #          ),
                    #          tags$hr(),
                    #          fluidRow(
                    #            column(width = 5, textOutput("text1"), tableOutput("heat.click")),
                    #            column(width = 5, textOutput("text2"), tableOutput("heat2.click"))
                    #          ))

                    # tabPanel("Differential Analysis",
                    #          uiOutput("differentialsText"),
                    #          tags$p(),
                    #          fluidRow(
                    #            column(width = 3, selectizeInput("differentialsCluster1Sel",
                    #                                             label    = "Cluster 1",
                    #                                             choices  = NULL,
                    #                                             multiple = FALSE)),
                    #            column(width = 3, selectizeInput("differentialsCluster2Sel",
                    #                                             label    = "Cluster 2",
                    #                                             choices  = NULL,
                    #                                             multiple = FALSE)),
                    #            column(width = 6,
                    #                   tags$label(HTML("&zwnj;")), tags$br(),
                    #                   create_bs_button(id    = "diffCalculateBtn",
                    #                                    label = "Calculate Differentials",
                    #                                    width = "60%",
                    #                                    disabled = TRUE))),
                    #          tags$hr(),
                    #          tags$h4(align = "center", uiOutput("differentialsTableTitle")),
                    #          uiOutput("differentialsTableAlternativeText"),
                    #          div(heatmap_downloadableTableUI("differentialsTable", list("csv", "tsv"),
                    #                                          "Download table data",
                    #                                          contentHeight = "200px"),
                    #              style = "padding-left:10%;padding-right:10%;"))
                     )
             )

# -- Register Elements in the ORDER SHOWN in the UI
add_ui_body(list(loading_modal, file_error_modal, body1))
