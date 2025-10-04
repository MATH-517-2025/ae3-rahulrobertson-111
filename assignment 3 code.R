library(shiny)
library(ggplot2)
library(viridis)
library(dplyr)

### ---- Core function ----
run_blockwise_analysis_2 <- function(n, N, alpha, beta, sd_noise = 1) {
  X <- rbeta(n, alpha, beta)
  m <- function(x) sin(1 / (x / 3 + 0.1))
  Y <- m(X) + rnorm(n, mean = 0, sd = sd_noise)
  
  m_double_prime <- function(x) {
    inner <- x / 3 + 0.1
    (2/9)*inner^(-3)*cos(1/inner) - (1/9)*inner^(-4)*sin(1/inner)
  }
  
  block_size <- floor(n / N)
  block_id <- rep(1:N, each = block_size)
  if(length(block_id) < n) block_id <- c(block_id, rep(N, n - length(block_id)))
  block_id <- block_id[1:n]
  
  XY <- data.frame(X = X, Y = Y)
  blocks <- split(XY, block_id)
  
  fits <- lapply(blocks, function(df) {
    if(nrow(df) >= 5) lm(Y ~ poly(X, 4, raw=TRUE), data=df) else NULL
  })
  valid_blocks <- which(!sapply(fits, is.null))
  fits <- fits[valid_blocks]
  coefs_list <- lapply(fits, coef)
  
  block_map <- rep(NA, N)
  block_map[valid_blocks] <- seq_along(valid_blocks)
  block_id_mapped <- block_map[block_id]
  
  m_hat_functions <- lapply(coefs_list, function(b) {
    function(x) b[1] + b[2]*x + b[3]*x^2 + b[4]*x^3 + b[5]*x^4
  })
  m_hat_primeprime_functions <- lapply(coefs_list, function(b) {
    function(x) 2*b[3] + 6*b[4]*x + 12*b[5]*x^2
  })
  
  theta_hat <- mean(sapply(1:n, function(i) {
    j <- block_id_mapped[i]
    if(!is.na(j)) m_hat_primeprime_functions[[j]](X[i])^2 else NA
  }), na.rm=TRUE)
  
  sigma_hat <- mean(sapply(1:n, function(i) {
    j <- block_id_mapped[i]
    if(!is.na(j)) (Y[i] - m_hat_functions[[j]](X[i]))^2 else NA
  }), na.rm=TRUE)
  
  h_hat_amise <- n^(-1/5) * ((35 * sigma_hat)/theta_hat)^(1/5)
  theta_true <- integrate(function(x) (m_double_prime(x))^2 * dbeta(x, alpha, beta), 0, 1)$value
  h_amise <- n^(-1/5) * (35 / theta_true)^(1/5)
  
  RSS <- sum(sapply(1:n, function(i) {
    j <- block_id_mapped[i]
    if(!is.na(j)) (Y[i] - m_hat_functions[[j]](X[i]))^2 else 0
  }))
  
  list(
    h_hat_amise=h_hat_amise,
    h_amise=h_amise,
    theta_hat=theta_hat,
    sigma_hat=sigma_hat,
    X=X,
    Y=Y,
    RSS=RSS,
    block_id=block_id_mapped,
    N=N,
    n=n
  )
}

choose_optimal_N <- function(n, alpha, beta, sd_noise = 1, N_candidates = 1:50) {
  # Compute N_max based on Ruppert et al. (1995)
  N_max <- max(min(floor(n / 20), 5), 1)
  
  # Compute RSS for N_max first
  res_max <- run_blockwise_analysis_2(n, N_max, alpha, beta, sd_noise)
  RSS_Nmax <- res_max$RSS
  denom <- RSS_Nmax / (n - 5 * N_max)
  
  # Compute criterion for all candidate N
  Cp_vals <- sapply(N_candidates, function(Nc) {
    res <- run_blockwise_analysis_2(n, Nc, alpha, beta, sd_noise)
    RSS_N <- res$RSS
    Cp <- RSS_N / denom - (n - 10 * Nc)
    return(Cp)
  })
  
  # Find optimal N
  N_opt <- N_candidates[which.min(Cp_vals)]
  
  list(N_candidates = N_candidates, Cp_vals = Cp_vals, N_opt = N_opt)
}
### ---- Shiny App ----
ui <- fluidPage(
  titlePanel("AMISE Bandwidth Explorer"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n","Sample size n:",min=100,max=10000,value=1000,step=100,animate=TRUE),
      numericInput("N","Number of blocks N:",value=10,min=1,step=1),
      sliderInput("alpha","Beta α:",min=0.1,max=10,value=1,step=0.1,animate=TRUE),
      sliderInput("beta","Beta β:",min=0.1,max=10,value=2,step=0.1,animate=TRUE),
      sliderInput("sd_noise","Noise SD:",min=0,max=5,value=1,step=0.1,animate=TRUE)
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Single run", plotOutput("h_plot"), verbatimTextOutput("summary_text")),
        tabPanel("error vs n", plotOutput("hn_plot")),
        tabPanel("Error vs N", plotOutput("errN_plot")),
        tabPanel("ĥ vs N", plotOutput("hhatN_plot")),
        tabPanel("ĥ vs n", plotOutput("hhat_vs_n_plot")),
        tabPanel("Sample distribution", plotOutput("sample_plot")),
        tabPanel("True vs Estimated h", plotOutput("h_vs_hhat_plot")),
        tabPanel("2D Histogram", plotOutput("hist2d_plot")),
        tabPanel("θ̂ & σ̂ vs N", plotOutput("theta_sigma_vs_N_plot")),
        tabPanel("Mallow’s Cp", plotOutput("cp_plot"))
      )
    )
  )
)

server <- function(input, output) {
  
  results <- reactive({
    run_blockwise_analysis_2(input$n, input$N, input$alpha, input$beta, input$sd_noise)
  })
  
  output$h_plot <- renderPlot({
    res <- results()
    df <- data.frame(Type=c("h_amise (true)","h_hat_amise"),Value=c(res$h_amise,res$h_hat_amise))
    ggplot(df,aes(x=Type,y=Value,fill=Type)) + geom_bar(stat="identity") +
      scale_fill_manual(values=c("blue","red")) + labs(title=paste("True vs estimated AMISE bandwidth (N =",res$N,")")) +
      theme_minimal()
  })
  
  output$summary_text <- renderPrint({
    res <- results()
    cat("True h_amise:",res$h_amise,"\n")
    cat("Estimated h_hat_amise:",res$h_hat_amise,"\n")
    cat("Theta_hat:",res$theta_hat,"\n")
    cat("Sigma_hat:",res$sigma_hat,"\n")
    cat("Chosen N:",res$N,"(manual)\n")
  })
  
  output$hn_plot <- renderPlot({
    # Create a fine grid of n values
    n_vals <- seq(100, 10000, by = 100)
    
    # Compute squared error for each n
    sq_errors <- sapply(n_vals, function(nn) {
      res <- run_blockwise_analysis_2(nn, input$N, input$alpha, input$beta, input$sd_noise)
      (res$h_amise - res$h_hat_amise)^2
    })
    
    df <- data.frame(n = n_vals, error = sq_errors)
    
    ggplot(df, aes(x = n, y = error)) +
      geom_line(color = "blue") +
      geom_point(data = df[seq(1, nrow(df), by = 10), ], aes(x = n, y = error), color = "red") + # optional sparser points
      labs(title = "Squared Error vs n", x = "n", y = "(h - h_hat)^2") +
      theme_minimal()
  })
  
  output$errN_plot <- renderPlot({
    N_vals <- c(1,5,10,20,50)
    sq_errors <- sapply(N_vals,function(Nc){
      res <- run_blockwise_analysis_2(input$n, Nc, input$alpha, input$beta, input$sd_noise)
      (res$h_amise - res$h_hat_amise)^2
    })
    df <- data.frame(N=N_vals,error=sq_errors)
    ggplot(df,aes(x=N,y=error)) + geom_line(color="blue") + geom_point(color="red") +
      labs(title="Squared Error vs N",x="N",y="(h - h_hat)^2") + theme_minimal()
  })
  
  output$hhatN_plot <- renderPlot({
    N_vals <- c(1,5,10,20,50)
    hhat <- sapply(N_vals,function(Nc){
      run_blockwise_analysis_2(input$n,Nc,input$alpha,input$beta,input$sd_noise)$h_hat_amise
    })
    df <- data.frame(N=N_vals,h_hat=hhat)
    ggplot(df,aes(x=N,y=h_hat)) + geom_line(color="red") + geom_point(color="red") +
      labs(title="h_hat vs N",x="N",y="h_hat") + theme_minimal()
  })
  
  output$hhat_vs_n_plot <- renderPlot({
    # Create a fine grid of n values
    n_vals <- seq(100, 10000, by = 100)
    
    # Compute estimated h_hat_amise for each n
    hhat <- sapply(n_vals, function(nn) {
      run_blockwise_analysis_2(nn, input$N, input$alpha, input$beta, input$sd_noise)$h_hat_amise
    })
    
    df <- data.frame(n = n_vals, h_hat = hhat)
    
    ggplot(df, aes(x = n, y = h_hat)) +
      geom_line(color = "red") +
      geom_point(data = df[seq(1, nrow(df), by = 10), ], aes(x = n, y = h_hat), color = "red") + # optional sparser points
      labs(title = "Estimated h_hat_amise vs n", x = "n", y = "h_hat_amise") +
      theme_minimal()
  })
  
  output$sample_plot <- renderPlot({
    res <- results()
    ggplot(data.frame(X=res$X),aes(x=X)) + geom_histogram(bins=30,fill="skyblue",color="black") +
      labs(title="Distribution of X ~ Beta(α,β)",x="X",y="Count") + theme_minimal()
  })
  
  output$h_vs_hhat_plot <- renderPlot({
    N_vals <- c(1,5,10,20,50)
    results_list <- lapply(N_vals,function(Nc) run_blockwise_analysis_2(input$n,Nc,input$alpha,input$beta,input$sd_noise))
    df <- data.frame(N=N_vals,
                     h_amise=sapply(results_list,function(r) r$h_amise),
                     h_hat_amise=sapply(results_list,function(r) r$h_hat_amise))
    ggplot(df) + 
      geom_line(aes(x=N,y=h_amise,color="h_amise"),size=1) +
      geom_point(aes(x=N,y=h_amise,color="h_amise")) +
      geom_line(aes(x=N,y=h_hat_amise,color="h_hat_amise"),size=1,linetype="dashed") +
      geom_point(aes(x=N,y=h_hat_amise,color="h_hat_amise")) +
      scale_color_manual(values=c("h_amise"="blue","h_hat_amise"="red")) +
      labs(title="True vs Estimated h",x="N",y="Bandwidth") + theme_minimal()
  })
  
  output$hist2d_plot <- renderPlot({
    res <- results()
    ggplot(data.frame(X=res$X,Y=res$Y),aes(x=X,y=Y)) + geom_bin2d(bins=30) +
      scale_fill_viridis() + labs(title="2D Histogram of (X,Y)") + theme_minimal()
  })
  
  output$theta_sigma_vs_N_plot <- renderPlot({
    N_vals <- 1:50
    results_list <- lapply(N_vals,function(Nc) run_blockwise_analysis_2(input$n,Nc,input$alpha,input$beta,input$sd_noise))
    df <- data.frame(N=N_vals,
                     theta_hat=sapply(results_list,function(r) r$theta_hat),
                     sigma_hat=sapply(results_list,function(r) r$sigma_hat))
    ggplot(df,aes(x=N)) +
      geom_line(aes(y=theta_hat,color="theta_hat"),size=1) +
      geom_point(aes(y=theta_hat,color="theta_hat")) +
      geom_line(aes(y=sigma_hat,color="sigma_hat"),size=1,linetype="dashed") +
      geom_point(aes(y=sigma_hat,color="sigma_hat")) +
      scale_color_manual(name="Estimate",values=c("theta_hat"="blue","sigma_hat"="red")) +
      labs(title=paste("θ̂ & σ̂ vs N (n =",input$n,")"),x="N",y="Value") + theme_minimal()
  })
  output$cp_plot <- renderPlot({
    opt <- choose_optimal_N(input$n,input$alpha,input$beta,input$sd_noise)
    df <- data.frame(N=opt$N_candidates,F=opt$Cp_vals)
    ggplot(df,aes(x=N,y=F)) +
      geom_line(color="purple") + geom_point(color="black") +
      geom_vline(xintercept=opt$N_opt, linetype="dashed", color="red") +
      labs(title=paste("Mallow's Cp criterion (Optimal N =", opt$N_opt,")"),
           x="N",y="F(N)")
  })
}

shinyApp(ui=ui,server=server)

