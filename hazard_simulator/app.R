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

CHAINS <- 1
CORES <- 2
ITER <- 500
SEED <- 42

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Survival Model Explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput("hazard_rate",
                        "Constant hazard rate (per 1000 units time): ",
                        min = 1,
                        max = 50,
                        value = 20
            ),
            shiny::actionButton('Add model', inputId = 'submit'),
            actionButton('Clear all', inputId = 'clear')
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
    fit_hazard_rates <- reactiveValues(data = NULL)
    
    # clear models 
    observeEvent(input$clear, {
        prior_surv$data <- NULL
        prior_haz$data <- NULL
        prior_cumhaz$data <- NULL
        fit_hazard_rates$data <- NULL
    })
    
    # add prior predictive
    observeEvent(input$submit, {
        hazard_rate <- isolate(input$hazard_rate)
        if (!hazard_rate %in% isolate(fit_hazard_rates$data)) {
            prior_model <- rstanarm::stan_surv(
                formula = Surv(time = end, event = event) ~ 1, # fit an intercept only model
                data = events,
                basehaz = "exp",
                prior_PD = TRUE,
                chains = CHAINS,
                cores = CORES,
                iter = ITER,
                seed = SEED,
                prior_intercept = normal(log(hazard_rate/1000), 0.0001) # prior on the log hazard
                )
            
            model_description <- glue::glue('Hazard rate = {scales::comma(hazard_rate)}/1000')
            # update prior predicted values
            prior_cumhaz$data <- posterior_survfit(prior_model,
                                                   type = 'cumhaz', 
                                                   newdata = events, prob = 0.5) %>%
                dplyr::group_by(time) %>%
                dplyr::summarise(median = median(median),
                                 ci_lb = median(ci_lb),
                                 ci_ub = median(ci_ub)) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(model = model_description) %>% 
                dplyr::bind_rows(isolate(prior_cumhaz$data))
            
            prior_surv$data <- posterior_survfit(prior_model,
                                                 type = 'surv', 
                                                 newdata = events, prob = 0.5) %>%
                dplyr::group_by(time) %>%
                dplyr::summarise(median = median(median),
                                 ci_lb = median(ci_lb),
                                 ci_ub = median(ci_ub)) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(model = model_description) %>%
                dplyr::bind_rows(isolate(prior_surv$data))
            
            prior_haz$data <- posterior_survfit(prior_model,
                                                type = 'haz', 
                                                newdata = events, prob = 0.5) %>%
                dplyr::group_by(time) %>%
                dplyr::summarise(median = median(median),
                                 ci_lb = median(ci_lb),
                                 ci_ub = median(ci_ub)) %>%
                dplyr::ungroup() %>%
                dplyr::mutate(model = model_description) %>%
                dplyr::bind_rows(isolate(prior_haz$data))
        }
        fit_hazard_rates$data <- c(isolate(fit_hazard_rates$data), hazard_rate)
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
                geom_lineribbon(alpha = 0.4) + 
                #theme(legend.position = 'none') +
                scale_y_continuous('Hazard rate (per 1000 units time)') +
                ggtitle('Prior predicted hazard over time') +
                labs(caption = stringr::str_c('Showing median + 50% CrI',
                                              sep = '\n'))
        }
    })
    
    # plot cumulative hazard
    output$cumhazardPlot <- renderPlot({
        if (!is.null(prior_cumhaz$data) && nrow(prior_cumhaz$data) > 0 ) {
            prior_cumhaz$data %>% 
                ggplot(.,
                       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = model, colour = model, fill = model)) + 
                geom_lineribbon(alpha = 0.4) + 
                scale_y_continuous('Cumulative hazard') +
                ggtitle('Prior predicted cumulative hazard over time') +
                labs(caption = stringr::str_c('Showing median + 50% CrI',
                                              sep = '\n'))
        }
    })
    
    # plot survival
    output$survivalPlot <- renderPlot({
        if (!is.null(prior_surv$data) && nrow(prior_surv$data) > 0) {
            prior_surv$data %>% 
                ggplot(.,
                       aes(x = time, y = median, ymin = ci_lb, ymax = ci_ub, group = model, colour = model, fill = model)) + 
                geom_lineribbon(alpha = 0.4) + 
                scale_y_continuous('Predicted survival (%)', labels = scales::percent) +
                ggtitle('Prior predicted survival over time') +
                labs(caption = stringr::str_c('Showing median + 50% CrI',
                                              sep = '\n'))
        }
    })    
}

# Run the application 
shinyApp(ui = ui, server = server)
