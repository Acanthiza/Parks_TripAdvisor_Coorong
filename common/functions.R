
  # glm 'chi-square'

  success_failure_GLM <- function(df){
    
    res <- list()
    
    var1 <- names(df)[1]
    var2 <- names(df)[2]
    
    dfSetup <- df %>%
      dplyr::group_by(!!ensym(var2)) %>%
      dplyr::mutate(levels = n()) %>%
      dplyr::ungroup() %>%
      dplyr::filter(levels >= 2/3*max(levels)) %>%
      dplyr::mutate(!!var2 := factor(!!ensym(var2))
                    , trials = success + failure
                    ) %>%
      dplyr::group_by(!!ensym(var1),trials) %>%
      tidyr::expand(!!ensym(var2)) %>%
      dplyr::ungroup()
    
    outcome <- df %>%
      dplyr::right_join(dfSetup) %>%
      dplyr::mutate(success = ifelse(is.na(success),0,success)
                    , failure = ifelse(is.na(failure),trials,failure)
                    , !!var2 := factor(!!ensym(var2))
                    , var1 = factor(!!ensym(var1))
                    , var2 = factor(!!ensym(var2))
                    )
    
    res$mod <- stan_glm(cbind(success,trials-success) ~ var1*var2
                    , data = outcome
                    , family = binomial()
                    , iter = 5000
                    )
    
    res$modFit <- pp_check(res$mod)
    
    res$modRhat <- plot(res$mod, "rhat_hist")
    
    res$modTrace <- stan_trace(res$mod)
    
    res$mod2d <- pp_check(res$mod, plotfun = "stat_2d", stat = c("mean", "sd"))
    
    # Use the model to predict results over variables of interest
    res$modPred <- outcome %>%
      dplyr::group_by(var1,var2) %>%
      dplyr::summarise(success = 0, trials = 100) %>% # use 100 trials to give results as percentages
      dplyr::ungroup() %>%
      dplyr::mutate(col = row.names(.)) %>%
      dplyr::left_join(as_tibble(posterior_predict(res$mod
                                                   , newdata = .
                                                   )
                                 ) %>%
                         tibble::rownames_to_column(var = "row") %>%
                         tidyr::gather(col,value,2:ncol(.))
                       ) %>%
      dplyr::left_join(outcome %>%
                         dplyr::select(!!ensym(var1),!!ensym(var2),var1,var2) %>%
                         unique()
                       )
    
    # summarise the results
    res$modRes <- as_tibble(res$modPred) %>%
      dplyr::group_by(!!ensym(var1),!!ensym(var2)) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(res$mod))
                       , modMedian = quantile(value,0.5)
                       , modMean = mean(value)
                       , modci90lo = quantile(value, 0.025)
                       , modci90up = quantile(value, 0.975)
                       , ci = modci90up-modci90lo
                       , text = paste0(round(modMedian,2)," (",round(modci90lo,2)," to ",round(modci90up,2),")")
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(outcome) %>%
      dplyr::mutate_if(is.numeric,round,2)
    
    
    res$modPlotRidges <- ggplot(res$modPred, aes(value,!!ensym(var1),fill=!!ensym(var1))) +
      ggridges::geom_density_ridges(alpha = 0.5) +
      facet_wrap(~get("var2"),scales="free") +
      scale_fill_viridis_d()
    
    
    res$modDiffVar1 <- res$modPred %>%
      (function(x) x %>% dplyr::left_join(x %>% dplyr::select(var1,var2,row,value) %>% dplyr::rename(var1b = var1, value_2 = value))) %>%
      dplyr::filter(var1 != var1b) %>%
      dplyr::mutate(diff = value-value_2
                    , comparison = map2_chr(var1,var1b,~paste(sort(c(.x,.y))[1],sort(c(.x,.y))[2]))
                    ) %>%
      dplyr::group_by(comparison,var2,row) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    
    res$modDiffRes <- res$modDiffVar1 %>%
      dplyr::group_by(var1,var2,var1b) %>%
      dplyr::summarise(n = n()
                       , nCheck = nrow(as_tibble(res$mod))
                       #, value = median(value)
                       , value = 100*sum(diff > 0)/n
                       #, `l = r` = 100*sum(value ==0)/n
                       ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(alpha = 1.5*abs(value-50)
                    , colour = if_else(value > 50,"Greater than 50%","Less than 50%")
                    , text = paste0(round(value,0),"% chance that ",var1," reviewers use\n",var2," more than ",var1b," reviewers")
                    , var1Comparison = paste0(var1," vs ",var1b)
                    , var1Comparison = fct_relevel(var1Comparison, grep("SA",unique(var1Comparison),value=TRUE))
                    ) %>%
      dplyr::mutate_if(is.numeric,round,2) %>%
      dplyr::arrange(var2,var1,var1b)
    
    
    res$modPlotRidgesDiff <- ggplot(res$modDiffVar1, aes(diff,paste0(get("var1")," vs ",get("var1b")))) +
      ggridges::geom_density_ridges(alpha = 0.5) +
      geom_vline(aes(xintercept = 0)) +
      facet_wrap(~get("var2"),scales="free") +
      scale_fill_viridis_d() +
      labs(y = "Comparison"
           , x = "Difference"
           )
      
    
    res$modPlotDiff <- ggplot(res$modDiffRes,aes(var1Comparison,var2,fill=colour,alpha=alpha,label=text)) +
      geom_tile() +
      geom_text(size=2.5) +
      scale_fill_viridis_d() +
      scale_alpha_continuous(guide = FALSE
                             , range = c(0,0.5)
                             ) +
      scale_x_discrete(limits = rev(levels(var2))) +
      labs(subtitle = "50% represents equal chance. Further from 50% is less faded (greater difference)"
           , fill = "Likelihood"
           , x = "Comparison"
           , y = get("var2")
           )
    
    return(res)
    
  }

  # chi-square
  chi_square <- function(cont) {
    
    var1 <- names(cont)[1]
    var2 <- names(cont)[2]
    
    contingency <- cont %>%
      dplyr::mutate(var1 = factor(!!ensym(var1))
                    , var2 = factor(!!ensym(var2))
                    ) %>%
      tidyr::complete(var1,var2) %>%
      dplyr::mutate(!!var1 := if_else(is.na(!!ensym(var1)),var1,!!ensym(var1))
                    , !!var2 := if_else(is.na(!!ensym(var2)),var2,!!ensym(var2))
                    , var1 = fct_inorder(factor(as.character(var1)))
                    , var2 = fct_inorder(factor(as.character(var2)))
                    , var1No = as.factor(as.numeric(var1))
                    , var2No = as.factor(as.numeric(var2))
                    ) %>%
      replace(is.na(.), 0)
    
    chSq <- contingency %>%
      dplyr::select(var1No,var2No,n) %>%
      tidyr::spread(var2No,n,drop=TRUE) %>%
      as.data.frame %>%
      tibble::column_to_rownames(names(.)[1]) %>%
      chisq.test()
    
    chSqResidual <- chSq$residuals %>%
      data.frame %>%
      tibble::rownames_to_column("var1No") %>%
      tidyr::gather("var2No","residual",2:ncol(.)) %>%
      dplyr::mutate(var2No = gsub("X","",var2No))
    
    chSqVis <- data.frame(100*chSq$residuals^2/chSq$statistic) %>%
      data.frame %>%
      tibble::rownames_to_column("var1No") %>%
      tidyr::gather("var2No","contribution",2:ncol(.))%>%
      dplyr::mutate(var2No = gsub("X","",var2No)) %>%
      as_tibble() %>%
      dplyr::left_join(chSqResidual) %>%
      dplyr::left_join(contingency) %>%
      dplyr::mutate(per = 100*contribution/sum(contribution)
                    , text = paste0("n:",n,"\n",round(per,1),"%")
                    , direction = if_else(residual<0
                                          ,"less than expected"
                                          , if_else(residual>0
                                                    ,"more than expected"
                                                    ,"as expected"
                                          )
                    )
                    , label = paste0(var2, " occurs ", direction, " in ", var1)
      ) %>%
      dplyr::select(!!var1,!!var2,contribution,residual,n,per,text,label,direction)
    
    chSqPlot <- ggplot(chSqVis, aes(!!ensym(var1), fct_rev(!!ensym(var2)), fill = direction, alpha = contribution, label = text)) +
      geom_tile() +
      geom_text(size = 2) +
      guides(alpha = FALSE) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
      labs(subtitle = "Percentages are the percent contribution to overall chi-squared value"
           , y = var2
           , x = var1
      ) +
      scale_fill_viridis_d()
    
    chSqText <- paste0("(Chi-squared = ",round(chSq$statistic,1), ", df = ",chSq$parameter,", p <= ",round(chSq$p.value,4),")")
    
    doChSqPlot <- chSq$p.value<0.05
    
    chSqRes <- list(chSq=chSq,chSqVis=chSqVis,chSqPlot=chSqPlot,chSqText=chSqText,doChSqPlot=doChSqPlot)
    
  }  
  
  
  
# Create a colour palette for n groups

  col_pal <-  function(n) {
    if (n <= 8) {
      RColorBrewer::brewer.pal(n, "Set2")
    } else {
      hcl(h=seq(0,(n-1)/(n),length=n)*360,c=100,l=65,fixup=TRUE)
    }
  }
  

  # turn a vector into a comma separated list of values with a penultimate 'and'
  vec_to_sentence <- function(x,sep=",") {
    
    x[!is.na(x)] %>%
      paste(collapse = "JOINSRUS") %>%
      (function(x) if(sep == ";") {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", paste0(sep," and ")) %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
      } else {
        
        stringi::stri_replace_last_regex(x,"JOINSRUS", " and ") %>%
          str_replace_all("JOINSRUS",paste0(sep," "))
        
      }
      )
    
  }

# https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r

  numbers2words <- function(x){
    ## Function by John Fox found here:
    ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
    ## Tweaks by AJH to add commas and "and"
    helper <- function(x){

      digits <- rev(strsplit(as.character(x), "")[[1]])
      nDigits <- length(digits)
      if (nDigits == 1) as.vector(ones[digits])
      else if (nDigits == 2)
        if (x <= 19) as.vector(teens[digits[1]])
      else trim(paste(tens[digits[2]],
                      Recall(as.numeric(digits[1]))))
      else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and",
                                        Recall(makeNumber(digits[2:1]))))
      else {
        nSuffix <- ((nDigits + 2) %/% 3) - 1
        if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
        trim(paste(Recall(makeNumber(digits[
          nDigits:(3*nSuffix + 1)])),
          suffixes[nSuffix],"," ,
          Recall(makeNumber(digits[(3*nSuffix):1]))))
      }
    }
    trim <- function(text){
      #Tidy leading/trailing whitespace, space before comma
      text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
      #Clear any trailing " and"
      text=gsub(" and$","",text)
      #Clear any trailing comma
      gsub("\ *,$","",text)
    }
    makeNumber <- function(...) as.numeric(paste(..., collapse=""))
    #Disable scientific notation
    opts <- options(scipen=100)
    on.exit(options(opts))
    ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
              "eight", "nine")
    names(ones) <- 0:9
    teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
               "sixteen", " seventeen", "eighteen", "nineteen")
    names(teens) <- 0:9
    tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
              "ninety")
    names(tens) <- 2:9
    x <- round(x)
    suffixes <- c("thousand", "million", "billion", "trillion")
    if (length(x) > 1) return(trim(sapply(x, helper)))
    helper(x)
  }

# Fix html widget when not displayed
  widgetFix <- function(inputFile,outputFile){
    a = readLines(inputFile)
    output = paste(a,collapse="\n")
    output = gsub(">\n\n</div>","></div>",output)
    writeLines(output,outputFile)
    invisible(NULL)
  }

# Generate jenks breaks
# http://cainarchaeology.weebly.com/r-function-for-plotting-jenks-natural-breaks-classification.html
  
  plotJenks <- function(data, n=3, brks.cex=0.70, top.margin=10, dist=5){
    df <- data.frame(sorted.values=sort(data, decreasing=TRUE))
    Jclassif <- classIntervals(df$sorted.values, n, style = "jenks") #requires the 'classInt' package
    test <- jenks.tests(Jclassif) #requires the 'classInt' package
    df$class <- cut(df$sorted.values, unique(Jclassif$brks), labels=FALSE, include.lowest=TRUE) #the function unique() is used to remove non-unique breaks, should the latter be produced. This is done because the cut() function cannot break the values into classes if non-unique breaks are provided
    if(length(Jclassif$brks)!=length(unique(Jclassif$brks))){
      info <- ("The method has produced non-unique breaks, which have been removed. Please, check '...$classif$brks'")
    } else {info <- ("The method did not produce non-unique breaks.")}
    loop.res <- numeric(nrow(df))
    i <- 1
    repeat{
      i <- i+1
      loop.class <- classIntervals(df$sorted.values, i, style = "jenks")
      loop.test <- jenks.tests(loop.class)
      loop.res[i] <- loop.test[[2]]
      if(loop.res[i]>0.9999){
        break
      }
    }
    max.GoF.brks <- which.max(loop.res)
    plot(x=df$sorted.values, y=c(1:nrow(df)), type="b", main=paste0("Jenks natural breaks optimization; number of classes: ", n), sub=paste0("Goodness of Fit: ", round(test[[2]],4), ". Max GoF (", round(max(loop.res),4), ") with classes:", max.GoF.brks), ylim =c(0, nrow(df)+top.margin), cex=0.75, cex.main=0.95, cex.sub=0.7, ylab="observation index", xlab="value (increasing order)")
    abline(v=Jclassif$brks, lty=3, col="red")
    text(x=Jclassif$brks, y= max(nrow(df)) + dist, labels=sort(round(Jclassif$brks, 2)), cex=brks.cex, srt=90)
    results <- list("info"=info, "classif" = Jclassif, "breaks.max.GoF"=max.GoF.brks, "class.data" = df)
    return(results)
  }
  