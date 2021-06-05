# Aplicación 

# Load packages
library(shiny)
library(shinythemes)
library(lattice)
library(ggplot2)
library(readr)
library(kernlab)
library(caret)


# load model
load("LC50svmrnexlooc.rda")
load("LC50rvmrexcvclass.rda")

# Define UI 
ui <- fluidPage(
    theme = shinytheme("darkly"),
    # Application title
    titlePanel("Predicción y Clasificación del Valor CL50 de Daphnia Magna en moléculas orgánicas"),

    sidebarLayout(
        sidebarPanel(h4("Ingreso de datos:"),    
                     fluidRow(column(7, 
                                     br(), br(),
                                     textInput("molecula", label = h4("Nombre de la molécula"),value = "Ciclopentano"),
                                     br(), br(),
                                     numericInput(inputId = "MM", h4("Masa Molar"), value = 70.1, min = 0),
                                     br(), br(),
                                     numericInput(inputId = "TPSA.Tot", h4("TPSA.Tot"), value = 0, min = 0, max = 100),
                                     br(), br(),
                                     numericInput(inputId = "H.050", h4("H.050"), value = 0, min = 0, max = 100, step = 1),
                                     br(), br(),
                                     numericInput(inputId = "MLOGP", h4("MLOGP"), value = 2.746),
                                     br(), br(),
                                     numericInput(inputId = "RDCHI", h4("RDCHI"), value = 1.667, min = 0, max = 100),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     actionButton(inputId = "go", "Actualizar"),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br(),
                                     br(), br()
                    
        ))),
        
        mainPanel(
            br(), br(),
            h2("Test"),
            h4("Predicción del valor LC50 obtenido a través de pruebas de toxicidad aguda 
               prácticadas en neonatos de Daphnia Magna de menos de 24 horas de nacidos en pruebas de toxicidad 
               aguda de 48 horas."),
            h4("Moléculas con valores menores o iguales de 10 mg/L son consideradas tóxicas."),
            br(), 
            h4("Daphnia Magna", align = "left"),
            img(src = "DM.png", height = 183, width = 150),
            br(), br(), 
            strong("Para realizar la predicción de una molécula orgánica se deberán ingresar 
                   los siguientes descriptores moleculares:"),
            br(), br(),
            em("* TPSA.Tot: superficie polar topológica calculada mediante un método de 
               contribución que considera N, O, P y S."),
            br(), br(),
            em("* H.050: número de átomos de hidrógeno unidos a heteroátomos."),
            br(), br(),
            em("* MLOGP: coeficiente de partición octanol-agua."),
            br(), br(),
            em("* RDCHI: índice topológico que codifica información sobre el tamaño molecular 
               y la ramificación, sin considerar a los heteroátomos."),
            br(), br(),
            h2("Modelo"),
            p("Las predicciones se realizan a través del modelo máquina de vectores de soporte con kernel radial generado en r."),
            p("La base de datos usada para el entrenamiento y prueba fue tomada de:"),
            a(href = "https://michem.unimib.it/download/data/acute-aquatic-toxicity-to-daphnia-magna/", "michem.unimib.it", target="_blank"),
            br(), br(),
            hr(),

            mainPanel(
            h2("Validación"),
            verbatimTextOutput("validacion"),
            br()
            ), 

            mainPanel(
            h2("Resultados"),
            br(),
            h4('LC50 de:'), 
            verbatimTextOutput("molecula"),
            h4('en mol/L:'),
            verbatimTextOutput("LC50", placeholder = TRUE),
            h4('en mg/L:'), 
            verbatimTextOutput("LC50ppm", placeholder = TRUE),
            h4('Clasificación:'), 
            verbatimTextOutput("Toxicidad", placeholder = TRUE),
            br()
            )
            
        ))
)


# Define server logic 
server <- function(input, output) {

      
      output$molecula = renderText(input$molecula)

      validacion = eventReactive(input$go, {
     
           ifelse(input$TPSA.Tot < 0 | input$H.050 < 0 | input$RDCHI < 0 | input$H.050 != round(input$H.050), 
             "TPSA.Tot, H.050 y RDCHI siempre toman valores positivos o cero y H.050 siempre es entero",
               "Valores validados")
       
           })
       
       output$validacion = validacion
       
       LC50 = eventReactive(input$go, {
           if(validacion() != "Valores validados") 
              { print("Valores ingresados no válidos") 
                  } else { 
                   df = t(as.numeric(data.frame(input$TPSA.Tot, input$H.050, input$MLOGP, input$RDCHI)))
                   colnames(df) = c("TPSA.Tot", "H.050", "MLOGP", "RDCHI")
                   round(10^(-predict(modelo.svm.r.p.se, df)), 6)
                  }  
       
          })

       output$LC50 = LC50
       
       LC50ppm = eventReactive(input$go, { if(input$MM <= 0) 
       { print("La masa molar siempre es positiva y diferente de cero") 
         } else if (is.numeric(LC50()) == "FALSE"){print("Valores ingresados no válidos")
       } else { 
           1000*input$MM*LC50()
         }
       
          })
       
       output$LC50ppm = LC50ppm
       
       Toxicidad = eventReactive(input$go, { if(is.numeric(LC50ppm()) == "FALSE" ) 
       { print("No se han validado los datos o no ha ingresado la Masa Molar ") 
       } else { 
         df = t(as.numeric(data.frame(input$TPSA.Tot, input$H.050, input$MLOGP, input$RDCHI)))
         colnames(df) = c("TPSA.Tot", "H.050", "MLOGP", "RDCHI")
         predict(modelo.svm.r.p.es, df)
       }
         
       })
       
       output$Toxicidad = Toxicidad
     
}

# Run the application 
shinyApp(ui = ui, server = server)

