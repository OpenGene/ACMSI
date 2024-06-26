# 引入所需的包
library(shiny)
library(plotly)
library(shinydashboard)
library(shinycssloaders)

# 定义UI
ui <- dashboardPage(
  # 设置标题
  dashboardHeader(title = "4D-Parameters with AUC"),
  dashboardSidebar(),
  # 设置主页面
  dashboardBody(
    # 创建选项卡
    tabsetPanel(
      tabPanel("3D大图", 
      sliderInput("threshold", "AUC阈值", min = 0.7, max = 1, value = 0, step = 0.001),
      shinycssloaders::withSpinner(plotlyOutput("scatterPlot"))),
      tabPanel("4个子图",
        fluidRow(
          column(6, shinycssloaders::withSpinner(plotlyOutput("scatterPlot1"))),
          column(6, shinycssloaders::withSpinner(plotlyOutput("scatterPlot2"))),
          column(6, shinycssloaders::withSpinner(plotlyOutput("scatterPlot3"))),
          column(6, shinycssloaders::withSpinner(plotlyOutput("scatterPlot4")))
        )
      ),
      tabPanel("AUC趋势图",
        shinycssloaders::withSpinner(plotlyOutput("aucPlot"))),
      tabPanel("三维图",
        shinycssloaders::withSpinner(plotlyOutput("ddPlot")))
    )
  )
)


# 定义服务器逻辑
server <- function(input, output) {
  
  # 生成数据
  ym_space <- vector()
  epr_space <- vector()
  epm_space <- vector()
  npm_space <- vector()
  auc_history <- vector()

  # 运行Shiny应用
  con <- file("/haplox/users/wenjl/project/MSI/paras0525.txt", "r")
  while (length(line <- readLines(con, n = 1)) > 0) {
    values <- strsplit(line, ", ")[[1]]
    ylim <- as.integer(strsplit(values[1], ": ")[[1]][2])
    extra_peak_bott <- as.integer(strsplit(values[2], ": ")[[1]][2])
    extra_peak_morph <- as.integer(strsplit(values[3], ": ")[[1]][2])
    new_peak_morph <- as.integer(strsplit(strsplit(values[4], ": ")[[1]][2], " ")[[1]][1])
    auc <- as.numeric(strsplit(values[4], ": ")[[1]][3])

    # 将值添加到列表
    ym_space <- c(ym_space, ylim)
    epr_space <- c(epr_space, extra_peak_bott)
    epm_space <- c(epm_space, extra_peak_morph)
    npm_space <- c(npm_space, new_peak_morph)
    auc_history <- c(auc_history, auc)
  }
  close(con)
    # 根据阈值筛选数据

  ## normalise
  # scaled_epr_space <- (epr_space - min(epr_space)) / (max(epr_space) - min(epr_space)) * 100
  # scaled_epm_space <- (epm_space - min(epm_space)) / (max(epm_space) - min(epm_space)) * 100
  # scaled_npm_space <- (npm_space - min(npm_space)) / (max(npm_space) - min(npm_space)) * 100
  # scaled_ym_space <- (ym_space - min(ym_space)) / (max(ym_space) - min(ym_space)) * 100
  # log_scaled_auc_history <- (log(auc_history) - log(0.9)) / (log(1) - log(0.9))
  ## not normalise
  scaled_epr_space <- epr_space
  scaled_epm_space <- epm_space
  scaled_npm_space <- npm_space
  scaled_ym_space <- ym_space
  log_scaled_auc_history <- auc_history
  ## 点的大小
  normalized_ym_space <- 30*(scaled_ym_space - min(scaled_ym_space)) / (max(scaled_ym_space) - min(scaled_ym_space))

  filtered_data <- reactive({
    data <- data.frame(
      scaled_epr_space,
      scaled_epm_space,
      scaled_npm_space,
      scaled_ym_space,
      normalized_ym_space,
      log_scaled_auc_history
    )
    data[data$log_scaled_auc_history > input$threshold, ]
  })
  # 4D
  generateScatterPlot <- function(filtered_data) {
    plot_ly(
      data = filtered_data(),
      x = ~scaled_epr_space,
      y = ~scaled_epm_space,
      z = ~scaled_npm_space,
      color = ~log_scaled_auc_history,
      type = "scatter3d",
      mode = "markers"
    ) %>%
    add_markers(
      marker = list(size = ~normalized_ym_space, opacity = 1),
      text = ~paste("NPR: ", scaled_epr_space, "<br>",
                    "EPM: ", scaled_epm_space, "<br>",
                    "NPM: ", scaled_npm_space, "<br>",
                    "YM: ", scaled_ym_space),
      hoverinfo = "text",
      scene = "scene"
    ) %>%
      layout(
        title = "EPR-EPM-NPM (ym): 4D-Paras",
        scene = list(
          xaxis = list(title = "NPR", tickvals = list(), ticktext = list()),
          yaxis = list(title = "EPM", tickvals = list(), ticktext = list()),
          zaxis = list(title = "NPM", tickvals = list(), ticktext = list()),
          camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
          aspectratio = list(x = 1, y = 1, z = 1),
          margin = list(l = 0, r = 0, b = 0, t = 0)
        ),
        height = 800
      )
  }
  
  output$scatterPlot <- renderPlotly({
    # 动态阈值
    filtered_data <- reactive({
      data <- data.frame(
        scaled_epr_space,
        scaled_epm_space,
        scaled_npm_space,
        scaled_ym_space,
        normalized_ym_space,
        log_scaled_auc_history
      )
      data[data$log_scaled_auc_history > input$threshold, ]
    })
    
    generateScatterPlot(filtered_data)
  })
  
  # 3D-4小图
  sizes = 10
  scatterPlot1 <- plot_ly(
    x = scaled_epr_space,
    y = scaled_npm_space,
    z = scaled_ym_space,
    color = log_scaled_auc_history,
    marker = list(size = sizes),
    type = "scatter3d",
    mode = "markers"
  ) %>% 
    layout(scene = list(xaxis = list(title = "EPR"),
                        yaxis = list(title = "NPM"),
                        zaxis = list(title = "YM"),
                        camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
                        aspectratio = list(x = 1, y = 1, z = 1),
                        margin = list(l = 0, r = 0, b = 0, t = 0))) %>%
    ggplotly()

  # 创建散点图2
  scatterPlot2 <- plot_ly(
    x = scaled_epr_space,
    y = scaled_npm_space,
    z = scaled_epm_space,
    color = log_scaled_auc_history,
    marker = list(size = sizes),
    type = "scatter3d",
    mode = "markers"
  ) %>% 
    layout(scene = list(xaxis = list(title = "EPR"),
                        yaxis = list(title = "NPM"),
                        zaxis = list(title = "EPM"),
                        camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
                        aspectratio = list(x = 1, y = 1, z = 1),
                        margin = list(l = 0, r = 0, b = 0, t = 0))) %>%
    ggplotly()

  # 创建散点图3
  scatterPlot3 <- plot_ly(
    x = scaled_ym_space,
    y = scaled_epm_space,
    z = scaled_npm_space,
    color = log_scaled_auc_history,
    marker = list(size = sizes),
    type = "scatter3d",
    mode = "markers"
  ) %>% 
    layout(scene = list(xaxis = list(title = "YM"),
                        yaxis = list(title = "EPM"),
                        zaxis = list(title = "NPM"),
                        camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
                        aspectratio = list(x = 1, y = 1, z = 1),
                        margin = list(l = 0, r = 0, b = 0, t = 0))) %>%
    ggplotly()

  # 创建散点图4
  scatterPlot4 <- plot_ly(
    x = scaled_ym_space,
    y = scaled_epr_space,
    z = scaled_npm_space,
    color = log_scaled_auc_history,
    marker = list(size = sizes),
    type = "scatter3d",
    mode = "markers"
  ) %>% 
    layout(scene = list(xaxis = list(title = "YM"),
                        yaxis = list(title = "EPR"),
                        zaxis = list(title = "NPM"),
                        camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
                        aspectratio = list(x = 1, y = 1, z = 1),
                        margin = list(l = 0, r = 0, b = 0, t = 0))) %>%
    ggplotly()

  generateddPlot <- function(filtered_data) {
    plot_ly(
      data = filtered_data(),
      x = ~scaled_epr_space,
      y = ~scaled_epm_space,
      z = ~log_scaled_auc_history,
      type = "scatter3d",
      mode = "markers",
      marker = list(
        opacity = 0.5
      )
    ) %>%
    add_trace(
      type = "surface",
      z = ~log_scaled_auc_history,
      colorscale = "Viridis",
      cmin = min(filtered_data()$log_scaled_auc_history),
      cmax = max(filtered_data()$log_scaled_auc_history)
    ) %>%
    layout(
      title = "EPM-EPR",
      scene = list(
        xaxis = list(title = "EPR"),
        yaxis = list(title = "EPM"),
        zaxis = list(title = "AUC"),
        camera = list(eye = list(x = 1.87, y = 0.88, z = -0.64)),
        aspectratio = list(x = 1, y = 1, z = 1),
        margin = list(l = 0, r = 0, b = 0, t = 0)
      ),
      height = 800
    )
  }
    
  output$ddPlot <- renderPlotly({
    # 动态阈值
    filtered_data <- reactive({
      data <- data.frame(
        scaled_epr_space,
        scaled_epm_space,
        scaled_npm_space,
        scaled_ym_space,
        normalized_ym_space,
        log_scaled_auc_history
      )
      # data[data$log_scaled_auc_history > input$threshold, ]
    })
    
    generateddPlot(filtered_data)
  })

  # 输出3D图
  output$scatterPlot1 <- renderPlotly(scatterPlot1)
  output$scatterPlot2 <- renderPlotly(scatterPlot2)
  output$scatterPlot3 <- renderPlotly(scatterPlot3)
  output$scatterPlot4 <- renderPlotly(scatterPlot4)
  # 输出趋势图
  output$aucPlot <- renderPlotly({
    aucPlot <- plot_ly(
      x = seq_along(auc_history),
      y = auc_history,
      type = "scatter",
      mode = "lines",
      fill = "tozeroy",
      line = list(color = "blue")
    ) %>%
      layout(xaxis = list(title = "Iteration"),
            yaxis = list(title = "AUC"))
    
    return(aucPlot)
  })
}

shinyApp(ui = ui, server = server, options = list(host = "0.0.0.0", port = 9527))