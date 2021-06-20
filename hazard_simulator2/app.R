#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(rstanarm)
library(tidyverse)
library(tidybayes)
ggplot2::theme_set(tidybayes::theme_ggdist())
source(here::here('functions', 'helper_functions.R'))

CHAINS <- 1
CORES <- 2
ITER <- 500
SEED <- 42
PROB <- 0.5

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Survival Model Explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("basehaz",
                        "Baseline hazard type",
                        choices = c('exp', 'ms')),
            
            conditionalPanel(
                condition = "input.basehaz=='exp'",
                numericInput('intercept_mean', 'Hazard rate (per 1000)', value = '5'),
                textOutput('intercept_log_loc'),
                numericInput('intercept_sd', 'Hazard rate SD', value = '0.1')
            ),
            
            conditionalPanel(
                condition =  "input.basehaz=='ms'",
                numericInput('df', 'df', value = '0'),
                numericInput('intercept_mean', 'Intercept Mean', value = '0'), 
                numericInput('intercept_sd', 'Intercept SD', value = '1')
            ),
            
            shiny::actionButton('Add model', inputId = 'submit'),
            actionButton('Clear all', inputId = 'clear'),
            h5('Click "Add model" to plot predicted values')
        ),
        
        # Show plots of the simulated hazard, 
        mainPanel(
            tabsetPanel(
                tabPanel("hazard", 
                   plotOutput("hazardPlot")
                   ),
                tabPanel("cumulative hazard",
                         plotOutput("cumhazardPlot"),
                ),
                tabPanel("survival",
                         plotOutput("survivalPlot")
                ),
                tabPanel("RMST",
                         plotOutput("rmstPlot")
                )
            )
        ) # end mainPanel
    ) # end sidebarLayout
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    events <- tibble::tibble(id = seq_len(1), end = 350, event = FALSE)
    prior_surv <- reactiveValues(data = NULL)
    prior_haz <- reactiveValues(data = NULL) 
    prior_cumhaz <- reactiveValues(data = NULL)
    prior_rmst <- reactiveValues(data = NULL)

    # clear models 
    observeEvent(input$clear, {
        prior_surv$data <- NULL
        prior_haz$data <- NULL
        prior_cumhaz$data <- NULL
        prior_rmst$data <- NULL
    })
    
    # add prior predictive
    observeEvent(input$submit, {
        if (input$basehaz == 'exp') {
            hazard_rate <- isolate(input$intercept_mean)
            log_hazard_rate <- log(isolate(input$intercept_mean)/1000)
            log_hazard_sd <- isolate(input$intercept_sd)
            prior_model <- rstanarm::stan_surv(
                formula = Surv(time = end, event = event) ~ 1, # fit an intercept only model
                data = events,
                basehaz = "exp",
                prior_PD = TRUE,
                chains = CHAINS,
                cores = CORES,
                iter = ITER,
                seed = SEED,
                prior_intercept = normal(log_hazard_rate, log_hazard_sd) # prior on the log hazard
            )
            
            model_description <- glue::glue('Constant hazard at {scales::comma(hazard_rate)}/1000 (sd = {scales::comma(log_hazard_sd, accuracy = 0.1)})')
        }
        # update prior predicted values
        prior_cumhaz$data <- posterior_survfit(prior_model,
                                               type = 'cumhaz', 
                                               newdata = events, prob = PROB) %>%
            dplyr::group_by(time) %>%
            dplyr::summarise(median = median(median),
                             ci_lb = median(ci_lb),
                             ci_ub = median(ci_ub)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(model = model_description) %>% 
            dplyr::bind_rows(isolate(prior_cumhaz$data))
        
        prior_surv$data <- posterior_survfit(prior_model,
                                             type = 'surv', 
                                             newdata = events, prob = PROB) %>%
            dplyr::group_by(time) %>%
            dplyr::summarise(median = median(median),
                             ci_lb = median(ci_lb),
                             ci_ub = median(ci_ub)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(model = model_description) %>%
            dplyr::bind_rows(isolate(prior_surv$data))
        
        prior_haz$data <- posterior_survfit(prior_model,
                                            type = 'haz', 
                                            newdata = events, prob = PROB) %>%
            dplyr::group_by(time) %>%
            dplyr::summarise(median = median(median),
                             ci_lb = median(ci_lb),
                             ci_ub = median(ci_ub)) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(model = model_description) %>%
            dplyr::bind_rows(isolate(prior_haz$data))
        
        prior_rmst$data <- rmst_check(prior_model, tau = 300, newdata = events) %>%
            dplyr::mutate(model = model_description) %>%
            dplyr::bind_rows(isolate(prior_rmst$data))
    })

    # plot hazard
    output$hazardPlot <- renderPlot({
        if (!is.null(prior_haz$data) && nrow(prior_haz$data) > 0) {
            prior_haz$data %>% 
                dplyr::mutate(median = median*1000,
                              ci_lb = ci_lb*1000,
                              ci_ub = ci_ub*1000
                              ) %>%
                ggplot(.,
                       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = model, colour = model, fill = model)) + 
                geom_lineribbon(alpha = 0.2) + 
                #theme(legend.position = 'none') +
                scale_y_continuous('Hazard rate (per 1000 units time)') +
                ggtitle('Prior predicted hazard over time') +
                labs(caption = glue::glue('Showing median + {scales::percent(PROB)} CrI'))
        }
    })

    # plot RMST
    output$rmstPlot <- renderPlot({
        if (!is.null(prior_rmst$data) && nrow(prior_rmst$data) > 0) {
            ggplot(prior_rmst$data,
                   aes(x = `1`, y = model, colour = model, fill = model)) + 
                tidybayes::stat_dotsinterval() + 
                theme(legend.position = 'none') +
                scale_x_continuous('Restricted mean survival time (tau = 300)') +
                ggtitle('Prior RMST at t = 300')
        }
    })
    
    # plot cumulative hazard
    output$cumhazardPlot <- renderPlot({
        if (!is.null(prior_cumhaz$data) && nrow(prior_cumhaz$data) > 0 ) {
            prior_cumhaz$data %>% 
                ggplot(.,
                       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = model, colour = model, fill = model)) + 
                geom_lineribbon(alpha = 0.2) + 
                scale_y_continuous('Cumulative hazard') +
                ggtitle('Prior predicted cumulative hazard over time') +
                labs(caption = glue::glue('Showing median + {scales::percent(PROB)} CrI'))
        }
    })
    
    # plot survival
    output$survivalPlot <- renderPlot({
        if (!is.null(prior_surv$data) && nrow(prior_surv$data) > 0) {
            prior_surv$data %>% 
                ggplot(.,
                       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = model, colour = model, fill = model)) + 
                geom_lineribbon(alpha = 0.2) + 
                scale_y_continuous('Predicted survival (%)', labels = scales::percent) +
                ggtitle('Prior predicted survival over time') +
                labs(caption = glue::glue('Showing median + {scales::percent(PROB)} CrI'))
        }
    })    
}

# Run the application 
shinyApp(ui = ui, server = server)
