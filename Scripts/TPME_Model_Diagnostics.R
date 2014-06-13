
Generate_Model_Diagnsotics <- function(input_folder_path = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/", input_file = "Test",LS_Actor = 8, out_directory = "~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/",Thin_Itterations = 1, vocab = vocabulary,county_name = "Testing",skip_first = 0, Cluster_Integrated = T){
    #load current results
    #load("Current_Itteration_Results.Rdata")
    #for testing load("~/Dropbox/PINLab/Projects/Denny_Working_Directory/2011_Analysis_Output/Transylvania_Sample_10M_2011_3-13-14.Rdata")
    
    
    if(Cluster_Integrated){
        #plots to generate
        library(coda)
        #library(statnet)
        library(gregmisc)
        library(ggplot2)
        print("Loading Data...")
        load(paste(input_folder_path,input_file,".Rdata", sep = ""))
        
        print("Extracting Reduced Data")
        #extract current metropolis results
        #Return_List = Result
        #vocab = vocabulary
        #Thin_Itterations = 1
        #skip_first=2000
        first_return <- 13
        Topic_Model_Results <- Return_List[1:5]
        Model_Parameters <- Return_List[6:first_return]
        Cluster_Topic_Assignments <- Return_List[(first_return+1):(first_return+Model_Parameters[[2]])]
        Last_Cluster_Topic_Assignments <- unlist(Cluster_Topic_Assignments[Model_Parameters[[2]]])
        Metropolis_Results <- Return_List[(first_return+1+Model_Parameters[[2]]):length(Return_List)]
        #str(Cluster_Topic_Assignments)
        #free up memory
        #rm(Return_List)
        
        
#         to_return[0] = token_topic_assignment_list;
#         to_return[1] = topic_present_edge_counts;
#         to_return[2] = topic_absent_edge_counts;
#         to_return[3] = token_type_topic_counts;
#         to_return[4] = edge_topic_assignments;
#         to_return[5] = number_of_documents;
#         to_return[6] = number_of_outer_itterations;
#         to_return[7] = number_of_Gibbs_itterations;
#         to_return[8] = number_of_MH_itterations;
#         to_return[9] = number_of_clusters;
#         to_return[10] = cur_proposal_variances;
#         to_return[11] = cur_accept_rates;
#         to_return[12] = Topic_Model_Likelihoods;
            
            #plot(Model_Parameters[[8]],ylim = c(-6000, -3800),ylab = "Unnormalized Topic Model Log Likelihood", xlab= "Iteration", pch = 20, col = "blue")
        
        skip_first= skip_first+1
        #remove the first skip_first itterations of each sublist and recombine
        Itterations <- Model_Parameters[[4]]
        temp <- append(Metropolis_Results[skip_first:Itterations],
                       Metropolis_Results[(Itterations+skip_first):(2*Itterations)])
        temp <- append(temp,Metropolis_Results[(2*Itterations+skip_first):(3*Itterations)])
        temp <- append(temp,Metropolis_Results[(3*Itterations+skip_first):(4*Itterations)])
        temp <- append(temp,Metropolis_Results[(4*Itterations+skip_first):(5*Itterations)])
        temp <- append(temp,Metropolis_Results[(5*Itterations+skip_first):(6*Itterations)])
        Metropolis_Results <- temp

        Itterations <- Model_Parameters[[4]] - skip_first +1
        
        #thin out the data by taking every Thin_Itterations itteration for the metropolis step
        Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]
        
        #get model information and extract data
        Latent_Spaces <- length(Metropolis_Results[[(Itterations + 1)]][,1,1])
        Clusters <- length(Metropolis_Results[[1]])
        #Clusters <- 2
        Topics <- length(Cluster_Topic_Assignments[[1]])
        Actors <- length(Metropolis_Results[[(Itterations + 1)]][1,1,])
        Token_Topic_Assignments <- Topic_Model_Results[[1]]
        Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
        Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
        Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
        Edge_Topic_Assignments <- Topic_Model_Results[[5]]
        Proposal_Variances <- Model_Parameters[[6]][Model_Parameters[[2]],]
        Cluster_Topic_Assigns <- Cluster_Topic_Assignments[[Model_Parameters[[2]]]]
        
        #get the total number of tokens assigned to each topic
        Topic_Token_Totals <- apply(Word_Type_Topic_Counts,2,sum)
        #get the total number of present edges assigned to each topic
        Email_Assignments <- apply(Topic_Present_Edge_Counts,3,sum)
        #number of words
        num_words <- length(vocab[,1])
        
        #this list will be used to store lists of three vectors: word indicies in descending order of likelihood, word probability in topic, actual word. 
        Topic_Top_Words <- list()
        Top_Ten_Words <- rep("",Topics)
        Top_Four_Words <- rep("",Topics)
        
        
        for(i in 1:Topics){
            indicies <- order(Word_Type_Topic_Counts[,i],decreasing = TRUE)
            probabilities <- Word_Type_Topic_Counts[indicies,i]/Topic_Token_Totals[i]
            #reduce to only words with non-zero probability
            indicies <- indicies[which(probabilities > 0)]
            probabilities <- probabilities[which(probabilities > 0)]
            top_words <- vocab[indicies,]
            Topic_List <- list(top_words,probabilities,indicies)
            Topic_Top_Words <- append(Topic_Top_Words,list(Topic_List))
            #add top words
            words <- ""
            words4 <- ""
            for(j in 1:10){
                words <- paste(words,top_words[j],", ",sep = "")
            }
            for(j in 1:4){
                words4 <- paste(words4,top_words[j],", ",sep = "")
            }
            Top_Ten_Words[i] <- words
            Top_Four_Words[i] <- words4
        }
        data <- list(Top_Four_Words,Cluster_Topic_Assigns)
        save(data,file=paste(out_directory ,county_name,"Topic_Top_Words.Rdata", sep = ""))
        #print out the top words for each topic in each cluster
        for(i in 1:Clusters){
            print(i)
            inds <- which(data[[2]] == i)
            print(data[[1]][inds])
        }




        #this list will be used to store lists of three vectors: word indicies in descending order of likelihood, word probability in topic, actual word. 
        Cluster_Top_Words <- list()
        Cluster_Ten_Words <- rep("",Clusters)
        Cluster_Four_Words <- rep("",Clusters)
        
        for(k in 1:Clusters){
            
            probabilities <- rep(0,length(Word_Type_Topic_Counts[,1]))
            for(i in 1:Topics){
                if(Cluster_Topic_Assigns[i] == k){
                    probabilities <- probabilities + Word_Type_Topic_Counts[,i]
                }
            }
            #reduce to only words with non-zero probability
            indicies <- order(probabilities,decreasing = TRUE)
            indicies <- indicies[which(probabilities > 0)]
            probabilities <- probabilities[which(probabilities > 0)]
            top_words <- vocab[indicies,]
            Cluster_List <- list(top_words,probabilities,indicies)
            Cluster_Top_Words <- append(Cluster_Top_Words,list(Cluster_List))
            #add top words
            words <- ""
            words4 <- ""
            for(j in 1:5){
                words <- paste(words,top_words[j],", ",sep = "")
            }
            words <- paste(words," \n",sep = "")
            for(j in 6:10){
                words <- paste(words,top_words[j],", ",sep = "")
            }
            for(j in 1:4){
                words4 <- paste(words4,top_words[j],", ",sep = "")
            }
            Cluster_Ten_Words[k] <- words
            Cluster_Four_Words[k] <- words4
        }
        
        
        
        
        
        
        
        #display current accept rate
        start <- (3*Itterations) + 1
        end <- 4*Itterations
        Accept_Rates <- Metropolis_Results[start:end]
        accepted <- matrix(0, ncol = Clusters,nrow = Itterations)
        for(j in 1:Itterations){
            accepted[j,] <- Accept_Rates[[j]]
        }

        for(i in 1:Clusters){
            print(paste("Current accept rate for cluster", i, "is:", sum(accepted[,i])/Itterations))
        }
        
        
        #======= plotting function definitions =======#
        #function to plot intercept over time
        plot_intercepts <- function(intercept){
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[i]][intercept]
            }
            plot(intercepts, main = paste("Topic:",intercept,"Geweke:",round(geweke.diag(intercepts)$z,2)),ylab= "Value",pch = 20)
        }
        #function to plot betas over time
        plot_betas <- function(Cluster){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][Cluster,]
            }
            matplot(betas, main = paste("Cluster:",Cluster,"Geweke","\n MM:",round(geweke.diag(betas[,1])$z,2) , "MF:",round(geweke.diag(betas[,2])$z,2) ,"\n FM:",round(geweke.diag(betas[,3])$z,2) ,"FF:",round(geweke.diag(betas[,4])$z,2)),ylab= "Value",pch = 20)
        }
        #function to plot latent spaces for a given actor
        plot_LS_Positions <- function(topic){
            LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
            for(i in 1:Itterations){
                LSP[i,] <- Metropolis_Results[[Itterations+i]][,topic,LS_Actor]
            }
            corel <- cor(LSP[,1], LSP[,2])
            matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3),"\n Geweke - LS1:",round(geweke.diag(LSP[,1])$z,2),"LS2:",round(geweke.diag(LSP[,2])$z,2)),ylab= "Value",pch = 20)
        }


        plot_Cluster_present_edge_Network <- function(Cluster){
            #get all topics assigned to cluster
            cluster_indexes <- which(Last_Cluster_Topic_Assignments == Cluster)
            #create matrix to add to
            slice <- matrix(0, Actors,Actors)
            if(length(cluster_indexes) > 0){
                for(i in cluster_indexes){
                    slice <- slice + Topic_Present_Edge_Counts[,,i]
                }
            }
            coordinates <- cbind(Metropolis_Results[[2*Itterations]][1,Cluster,],Metropolis_Results[[2*Itterations]][2,Cluster,])
            coordinates <- as.data.frame(coordinates)
            
            gend <- Author_Attributes$Gender
            colors <- rep("", length(gend))
            for(l in 1:length(gend)){
                if(gend[l] == "M"){
                    colors[l] <- "purple"
                }else{
                    colors[l] <- "orange"
                }
            }
            
            plot(coordinates,main = paste("Topic:",Cluster,"Number of Edges:",sum(slice),"\n", Cluster_Ten_Words[Cluster]),pch = 20, col = "red")
            #add in lines between actors
            for(i in 1:Actors){
                for(j in 1:Actors){
                    if(slice[i,j] >0){
                        lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                    }
                }
            }
            points(coordinates,col = colors, pch = 19, cex = 2)
            
        }


        plot_Cluster_absent_edge_Network <- function(Cluster){
            #get all topics assigned to cluster
            cluster_indexes <- which(Last_Cluster_Topic_Assignments == Cluster)
            #create matrix to add to
            slice <- matrix(0, Actors,Actors)
            if(length(cluster_indexes) > 0){
                for(i in cluster_indexes){
                    slice <- slice + Topic_Absent_Edge_Counts[,,i]
                }
            }
            coordinates <- cbind(Metropolis_Results[[2*Itterations]][1,Cluster,],Metropolis_Results[[2*Itterations]][2,Cluster,])
            coordinates <- as.data.frame(coordinates)
            plot(coordinates,main = paste("Topic:",Cluster,"Number of Edges:",sum(slice),"\n", Cluster_Four_Words[Cluster]),pch = 20, col = "red")
            #add in lines between actors
            for(i in 1:Actors){
                for(j in 1:Actors){
                    if(slice[i,j] >0){
                        lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                    }
                }
            }
            
        }
        
        #make plots will all top words for one per page output
        plot_Topic_Network_Full <- function(topic){
            slice <- Topic_Present_Edge_Counts[,,topic]
            coordinates <- cbind(Metropolis_Results[[Itterations]][1,topic,],Metropolis_Results[[Itterations]][2,topic,])
            coordinates <- as.data.frame(coordinates)
            plot(coordinates,main = paste("Topic:",topic,"Number of Edges:",sum(slice),"\n", Cluster_Ten_Words[topic]),pch = 20, col = "red")
            #add in lines between actors
            for(i in 1:Actors){
                for(j in 1:Actors){
                    if(slice[i,j] >0){
                        lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                    }
                }
            }
            
        }
        
        #function to plot log likelihoods over time
        plot_likelihoods <- function(Cluster){
            likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
            for(i in 1:Itterations){
                likelihoods[i,1] <- Metropolis_Results[[4*Itterations+i]][Cluster]
                likelihoods[i,2] <- Metropolis_Results[[5*Itterations+i]][Cluster]
            }
            matplot(likelihoods, main = paste("Log Likelihoods of Current and Proposed Positions over Time For Cluster", Cluster, "\n Geweke Statistic:", round(geweke.diag(likelihoods[,1])$z,2)),ylab= "Value",pch = 20)
        }
        
        #function to plot log ratios vs log uniform draws over time
        plot_ratio_lud <- function(Cluster){
            #likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
            accept <- rep(0,Itterations)
            likelihoods <- rep(0,Itterations)
            lud <- rep(0,Itterations)
            colors <- rep("",Itterations)
            for(i in 1:Itterations){
                likelihoods[i] <- Metropolis_Results[[4*Itterations+i]][Cluster] - Metropolis_Results[[5*Itterations+i]][Cluster] #- Metropolis_Results[[6*Itterations+i]]
                accept[i] <- accepted[i,Cluster]
                if(accept[i] == 1){
                    col <- "blue"
                }else{
                    col <- "red"
                }
                colors[i] <- col
                
                accepted2 <- accept*likelihoods
                accepted2 <- accepted2[which(accept == 1)]
                
            }
            #fit <- lm(accepted~1:10000)
            print(paste("Average Accepted Prob:",mean(accepted)))
            len <- length(accepted)
            print(paste("Average Accepted Prob last 10 percent:",mean(accepted[(len - len/10):len])))
            print(paste("Average Log Ratio:",mean(likelihoods)))
            print(paste("Average Log Uniform Draw:",mean(lud)))
            scatter.smooth(x = 1:length(accepted), y = accepted, main = "Log Likelihood Ratio of Accepted Proposals",ylab= "Value",pch = 20)
            plot(likelihoods, main = "Red represents rejected proposals, blue represents accepted proposals",ylab= "Value",pch = 20, col = colors)
        }
        
        #function to plot mean beta values over time
        plot_beta_estimates <- function(Cluster){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][Cluster,]
            }
            mean_se <- as.data.frame(matrix(0,nrow=4, ncol = 2))
            mean_se <- cbind(c("M-M", "M-F","F-M", "F-F"),mean_se)
            names(mean_se) <- c("Tie","Parameter_Estimate","SE")
            for(j in 1:length(mean_se[,1])){
                mean_se[j,2] <- mean(betas[,j])
                mean_se[j,3] <- sd(betas[,j]) 
            }
            plot <- ggplot(mean_se, aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue"),colour="black",pch=c(21,22,23,24), size = 4)+ labs(title = paste("Cluster:",Cluster))
            return(list(plot))
        }
        # ================================================ #
        
        
        
        # ======= Generate Plots and save as PDFs =========#      
        make_plots <- function(Width, Height, parrow,parcol){
            #generate pdf of intercepts
            print("Plotting Intercepts...")
            pdf(file=paste(out_directory,county_name,"_Intercepts.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_intercepts)
            dev.off()
            
            #generate pdf of intercepts 
            print("Plotting Betas...")
            pdf(file=paste(out_directory,county_name,"_Betas.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_betas)
            dev.off()
            
            #generate pdf of intercepts 
            print("Plotting Latent Space Positions...")
            pdf(file=paste(out_directory,county_name,"_LS_Positions.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_LS_Positions)
            dev.off()
            
#             #generate network plots
#             print("Plotting Networks...")
#             pdf(file=paste(out_directory,county_name,"_Network_Plots.pdf",sep = ""), width = Width, height = Height)
#             par(mfrow= c(parrow,parcol))
#             sapply(1:Clusters,plot_Topic_Network)
#             dev.off()
            
            pdf(paste(out_directory,county_name,"_Likelihoods.pdf", sep = ""),height=Height,width=Width,pointsize=7)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_likelihoods) 
            dev.off()
            
            #generate network plots one per page
            print("Plotting Networks...")
            pdf(file=paste(out_directory,county_name,"_Cluster_Absent_Edge_Networks.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_Cluster_absent_edge_Network)
            dev.off()

            pdf(file=paste(out_directory,county_name,"_Cluster_Present_Edge_Networks.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Clusters,plot_Cluster_present_edge_Network)
            dev.off()
            
            
            #generate network plots
            print("Plotting Beta Estimates...")
            pdf(file=paste(out_directory,county_name,"_Beta_Estimates.pdf",sep = ""), width = Width, height = 5)
            par(mfrow= c(parrow,parcol))
            plots <- sapply(1:Clusters,plot_beta_estimates)
            #for(j in 1:length(plots)){
            #    print(plots[[j]])
            #}
            multiplot(plotlist = plots, cols = 5)
            dev.off()
        }
        
        #stolen from the R cookbook
        multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
            require(grid)
            
            # Make a list from the ... arguments and plotlist
            plots <- c(list(...), plotlist)
            
            numPlots = length(plots)
            
            # If layout is NULL, then use 'cols' to determine layout
            if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols), byrow= T)
            }
            
            if (numPlots==1) {
                print(plots[[1]])
                
            } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                    
                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                }
            }
        }
        
        #now actually generate the plots based on the number of topics     
        if(Clusters <= 4){
            make_plots(12,12,2,2)
        }else if(Clusters <= 9 ){
            make_plots(12,12,3,3)
        }else if(Clusters <= 10){
            make_plots(15,9,2,5)
        }else if(Clusters <= 20){
            make_plots(25,25,4,5)
        }else{
            make_plots(50,25,6,10)   
        }
        
        par(mfrow= c(1,1))
        
        
        # ======= Now generate top words and topic model diagnostics ====== #
        
        print("Generating Topic Model Output")
        
        
        
        #Establish Ordering 
        #ordering <- order(Email_Assignments, decreasing = FALSE)
        #Email_Assignments <- Email_Assignments[ordering]
        #Top_Ten_Words <- Top_Ten_Words[ordering]
        
        
        
        pdf(paste(out_directory,county_name,"_Top_Words.pdf", sep = ""),height=12,width=10,pointsize=7)
        par(mfrow= c(1,1))
        bp <- barplot2(Email_Assignments, beside = TRUE, horiz = TRUE, col = "lightblue1", border= "lightblue1", ylab = "Topic:", xlab = "Number of Edges Assigned to Topic",main = paste("Top Words by Number of Words Assigned to Topic for ",county_name,sep = "")) 
        text(0,bp,Top_Ten_Words,cex=1,pos=4)# label on the bars themselves 
        dev.off()
        
        
        
        pdf(paste(out_directory,county_name,"_Log_Ratio_LUD.pdf", sep = ""),height=8,width=12,pointsize=7)
        par(mfrow= c(5,2))
        sapply(1:Clusters,plot_ratio_lud) 
        dev.off()
        
        print("Outputing Geweke statistics..")
        
        intercept_list <- list()
        beta_list <- list()
        LS_list <- list()
        t= 1
        for(t in 1:Clusters){
            print(paste("Current Topic: ", t))
            #Intercepts
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[i]][t]
            }
            int <- as.numeric(geweke.diag(intercepts)$z)
            intercept_list <- append(intercept_list,int)
            
            #Betas
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][t,]
            }
            bets <- rep(0,4)
            for(i in 1:4){
                bets[i] <- as.numeric(geweke.diag(betas[,i])$z)
            } 
            beta_list <- append(beta_list,bets)
            
            #Latent Space Positions
            lat <- matrix(0,nrow=Actors, ncol= Latent_Spaces)
            for(a in 1:Actors){
                LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
                for(i in 1:Itterations){
                    LSP[i,] <- Metropolis_Results[[Itterations+i]][,t,a]
                }
                #calcaulte statistic
                for(j in 1:Latent_Spaces){
                    lat[a,j] <- as.numeric(geweke.diag(LSP[,j])$z)
                }
            }
            LS_list <- append(LS_list,lat)
            
        }
        
        ints <- length(which(abs(unlist(intercept_list)) < 1.645))
        print(paste("Percent Intercepts Converged ($z < 1.645$):",ints/length(intercept_list), "Number:", ints))
        
        
        bets <- length(which(abs(unlist(beta_list)) < 1.645))
        print(paste("Percent Betas Converged ($z < 1.645$):",bets/length(beta_list), "Number:", bets))
        
        ls <- length(which(abs(unlist(LS_list)) < 1.645))
        print(paste("Percent LS Positions Converged ($z < 1.645$):",ls/length(LS_list), "Number:", ls))
        
        total <- ints + bets +ls
        print(paste("Percent All Parameters Converged ($z < 1.645$):",total/(length(beta_list)+ length(LS_list) + length(intercept_list)), "Number:", total))
        
        
        #generate asortativity - token dataset
        print("Generating topic- token- assortativity dataset")
        get_beta_estimates <- function(topic){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
            }
            mean_se <- as.data.frame(matrix(0,nrow=4, ncol = 2))
            mean_se <- cbind(c("M-M", "M-F","F-M", "F-F"),mean_se)
            names(mean_se) <- c("Tie","Parameter_Estimate","SE")
            for(j in 1:length(mean_se[,1])){
                mean_se[j,2] <- mean(betas[,j])
                mean_se[j,3] <- sd(betas[,j]) 
            }
            return(mean_se[,2:3])
        }
        
        beta_averages <- matrix(0, ncol = 4,nrow = Clusters)
        beta_se <- matrix(0, ncol = 4,nrow = Clusters)
        for(i in 1:Clusters){
            result <- get_beta_estimates(i)
            beta_averages[i,] <- result[,1]
            beta_se[i,] <- result[,2]
            
        }
        
        
        Clusters_intercepts <- matrix(0,ncol = 2,nrow= Topics)
        for(t in 1:Clusters){
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[i]][t]
            }
            Clusters_intercepts[t,1] <- mean(intercepts)
            Clusters_intercepts[t,2] <- sd(intercepts)
        }
        
        topic_average_distances <- rep(0,Clusters)
        for(t in 1:Clusters){
            print(paste("Calculating topic latent space average distance:",t))
            #slice <- Topic_Present_Edge_Counts[,,t]
            coordinates <- matrix(0,ncol = Latent_Spaces,nrow = Actors)
            for(i in 1:Actors){
                tem<- rep(0 ,Itterations)
                for(k in 1:Itterations){
                    tem[k] <- Metropolis_Results[[Itterations+ k]][1,t,i]
                }
                coordinates[i,1] <- mean(tem)
                tem2<- rep(0 ,Itterations)
                for(k in 1:Itterations){
                    tem2[k] <- Metropolis_Results[[Itterations +k]][2,t,i]
                }
                coordinates[i,2] <- mean(tem2)
            }
            #coordinates <- cbind(mean(Metropolis_Results[1:Itterations][1,t,]),mean(Metropolis_Results[1:Itterations][2,t,]))
            dist <- 0
            for(a in 1:Actors){
                for(b in 1:Actors){
                    if(a != b){
                        for(l in 1:Latent_Spaces){
                            dis <- (coordinates[a,1] - coordinates[b,2])^2
                        }
                        dist <- dist + sqrt(dis)
                    }
                }
            }
            dist <- dist/(Actors*(Actors -1))
            topic_average_distances[t] <- dist
        }
        
        
        
        #generate dataset for analysis
        transpose <- t(Word_Type_Topic_Counts)
        
        probs <- matrix(0, ncol = num_words, nrow = Topics)
        for(i in 1:Topics){
            probs[i,] <- Word_Type_Topic_Counts[,i]/Topic_Token_Totals[i] 
        }
        #now combine over topic in cluster
        probs2 <- matrix(0, ncol = num_words, nrow = Clusters)
        for(j in 1:Topics){
            clust <- Cluster_Topic_Assigns[j]
            probs2[clust,] <- probs2[clust,] + probs[j,]
        }
            

        topic_data <- cbind(Clusters_intercepts,topic_average_distances,beta_averages,beta_se,transpose,probs2)
        topic_data <- as.data.frame(topic_data)
        
        prob_vocab <- as.vector(vocab[,1])
        append_prob <- function(str){
            return(paste("PR-",str,sep = ""))
        }
        temp <- as.vector(sapply(prob_vocab,append_prob))
        vec <- c("Intercept","Intercept-SE","Average_LS_Distance","MM", "MF","FM", "FF","MM-SE", "MF-SE","FM-SE", "FF-SE",as.vector(vocab[,1]) ,temp)
        names(topic_data) <- vec
        save(topic_data,file=paste(out_directory ,county_name,"Topic_Assortativity_Data.Rdata", sep = ""))
        
        
        
        
        
        
        
        
        
        
        
        
    }else{
        
        
        ########################################################################
        ################################ NO CLUSTERING #########################
        ########################################################################
        #plots to generate
        library(coda)
        #library(statnet)
        library(gregmisc)
        library(ggplot2)
        print("Loading Data...")
        load(paste(input_folder_path,input_file,".Rdata", sep = ""))
        
        print("Extracting Reduced Data")
        #extract current metropolis results
        Metropolis_Results <- Return_List[[1]]
        Topic_Model_Results <- Return_List[[2]]
        #free up memory
        rm(Return_List)
        
        
        
        
        
        
        skip_first= skip_first+1
        #remove the first skip_first itterations of each sublist and recombine
        Itterations <- length(Metropolis_Results)/7
        temp <- append(Metropolis_Results[skip_first:Itterations],
                       Metropolis_Results[(Itterations+skip_first):(2*Itterations)])
        temp <- append(temp,Metropolis_Results[(2*Itterations+skip_first):(3*Itterations)])
        temp <- append(temp,Metropolis_Results[(3*Itterations+skip_first):(4*Itterations)])
        temp <- append(temp,Metropolis_Results[(4*Itterations+skip_first):(5*Itterations)])
        temp <- append(temp,Metropolis_Results[(5*Itterations+skip_first):(6*Itterations)])
        temp <- append(temp,Metropolis_Results[(6*Itterations+skip_first):(7*Itterations)])
        Metropolis_Results <- temp
        
        #thin out the data by taking every Thin_Itterations itteration for the metropolis step
        Metropolis_Results <- Metropolis_Results[seq(1, length(Metropolis_Results),Thin_Itterations)]
        
        #get model information and extract data
        Itterations <- length(Metropolis_Results)/7
        Latent_Spaces <- length(Metropolis_Results[[1]][,1,1])
        Topics <- length(Metropolis_Results[[1]][1,,1])
        Actors <- length(Metropolis_Results[[1]][1,1,])
        Token_Topic_Assignments <- Topic_Model_Results[[1]]
        Topic_Present_Edge_Counts <- Topic_Model_Results[[2]]
        Topic_Absent_Edge_Counts <- Topic_Model_Results[[3]]
        Word_Type_Topic_Counts <- Topic_Model_Results[[4]]
        Edge_Topic_Assignments <- Topic_Model_Results[[5]]
        
        
        #get the total number of tokens assigned to each topic
        Topic_Token_Totals <- apply(Word_Type_Topic_Counts,2,sum)
        #get the total number of present edges assigned to each topic
        Email_Assignments <- apply(Topic_Present_Edge_Counts,3,sum)
        #number of words
        num_words <- length(vocab[,1])
        
        #this list will be used to store lists of three vectors: word indicies in descending order of likelihood, word probability in topic, actual word. 
        Topic_Top_Words <- list()
        Top_Ten_Words <- rep("",Topics)
        Top_Four_Words <- rep("",Topics)
        
        
        for(i in 1:Topics){
            indicies <- order(Word_Type_Topic_Counts[,i],decreasing = TRUE)
            probabilities <- Word_Type_Topic_Counts[indicies,i]/Topic_Token_Totals[i]
            #reduce to only words with non-zero probability
            indicies <- indicies[which(probabilities > 0)]
            probabilities <- probabilities[which(probabilities > 0)]
            top_words <- vocab[indicies,]
            Topic_List <- list(top_words,probabilities,indicies)
            Topic_Top_Words <- append(Topic_Top_Words,list(Topic_List))
            #add top words
            words <- ""
            words4 <- ""
            for(j in 1:10){
                words <- paste(words,top_words[j],", ",sep = "")
            }
            for(j in 1:4){
                words4 <- paste(words4,top_words[j],", ",sep = "")
            }
            Top_Ten_Words[i] <- words
            Top_Four_Words[i] <- words4
        }
        
        
        
        
        
        
        
        
        #display current accept rate
        start <- (3*Itterations) + 1
        end <- 4*Itterations
        print(paste("Current_Accept_Rate:", sum(unlist(Metropolis_Results[start:end]))/Itterations))
        
        #======= plotting function definitions =======#
        #function to plot intercept over time
        plot_intercepts <- function(intercept){
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[Itterations+i]][intercept]
            }
            plot(intercepts, main = paste("Topic:",intercept,"Geweke:",round(geweke.diag(intercepts)$z,2)),ylab= "Value",pch = 20)
        }
        #function to plot betas over time
        plot_betas <- function(topic){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
            }
            matplot(betas, main = paste("Topic:",topic,"Geweke","\n MM:",round(geweke.diag(betas[,1])$z,2) , "MF:",round(geweke.diag(betas[,2])$z,2) ,"\n FM:",round(geweke.diag(betas[,3])$z,2) ,"FF:",round(geweke.diag(betas[,4])$z,2)),ylab= "Value",pch = 20)
        }
        #function to plot latent spaces for a given actor
        plot_LS_Positions <- function(topic){
            LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
            for(i in 1:Itterations){
                LSP[i,] <- Metropolis_Results[[i]][,topic,LS_Actor]
            }
            corel <- cor(LSP[,1], LSP[,2])
            matplot(LSP, main = paste("Topic:",topic, "\n LS Correlations:", round(corel,3),"\n Geweke - LS1:",round(geweke.diag(LSP[,1])$z,2),"LS2:",round(geweke.diag(LSP[,2])$z,2)),ylab= "Value",pch = 20)
        }
        plot_Topic_Network <- function(topic){
            slice <- Topic_Present_Edge_Counts[,,topic]
            coordinates <- cbind(Metropolis_Results[[Itterations]][1,topic,],Metropolis_Results[[Itterations]][2,topic,])
            coordinates <- as.data.frame(coordinates)
            plot(coordinates,main = paste("Topic:",topic,"Number of Edges:",sum(slice),"\n", Top_Four_Words[topic]),pch = 20, col = "red")
            #add in lines between actors
            for(i in 1:Actors){
                for(j in 1:Actors){
                    if(slice[i,j] >0){
                        lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                    }
                }
            }
            
        }
        
        #make plots will all top words for one per page output
        plot_Topic_Network_Full <- function(topic){
            slice <- Topic_Present_Edge_Counts[,,topic]
            coordinates <- cbind(Metropolis_Results[[Itterations]][1,topic,],Metropolis_Results[[Itterations]][2,topic,])
            coordinates <- as.data.frame(coordinates)
            plot(coordinates,main = paste("Topic:",topic,"Number of Edges:",sum(slice),"\n", Top_Ten_Words[topic]),pch = 20, col = "red")
            #add in lines between actors
            for(i in 1:Actors){
                for(j in 1:Actors){
                    if(slice[i,j] >0){
                        lines(c(coordinates[i,1],coordinates[j,1]) , c(coordinates[i,2],coordinates[j,2]), col = "black")
                    }
                }
            }
            
        }
        
        #function to plot log likelihoods over time
        plot_likelihoods <- function(Itterations){
            likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
            for(i in 1:Itterations){
                likelihoods[i,1] <- Metropolis_Results[[4*Itterations+i]]
                likelihoods[i,2] <- Metropolis_Results[[5*Itterations+i]]
            }
            matplot(likelihoods, main = paste("Log Likelihoods of Current and Proposed Positions over Time \n Geweke Statistic:", round(geweke.diag(likelihoods[,1])$z,2)),ylab= "Value",pch = 20)
        }
        
        #function to plot log ratios vs log uniform draws over time
        plot_ratio_lud <- function(Itterations){
            #likelihoods <- matrix(0,ncol = 2, nrow = Itterations)
            accept <- rep(0,Itterations)
            likelihoods <- rep(0,Itterations)
            lud <- rep(0,Itterations)
            colors <- rep("",Itterations)
            for(i in 1:Itterations){
                likelihoods[i] <- Metropolis_Results[[4*Itterations+i]] - Metropolis_Results[[5*Itterations+i]] #- Metropolis_Results[[6*Itterations+i]]
                accept[i] <- Metropolis_Results[[3*Itterations+i]]
                lud[i] <- Metropolis_Results[[6*Itterations+i]]
                if(Metropolis_Results[[3*Itterations+i]] == 1){
                    col <- "blue"
                }else{
                    col <- "red"
                }
                colors[i] <- col
                
                accepted <- accept*likelihoods
                accepted <- accepted[which(accept == 1)]
                
            }
            #fit <- lm(accepted~1:10000)
            print(paste("Average Accepted Prob:",mean(accepted)))
            len <- length(accepted)
            print(paste("Average Accepted Prob last 10 percent:",mean(accepted[(len - len/10):len])))
            print(paste("Average Log Ratio:",mean(likelihoods)))
            print(paste("Average Log Uniform Draw:",mean(lud)))
            scatter.smooth(x = 1:length(accepted), y = accepted, main = "Log Likelihood Ratio of Accepted Proposals",ylab= "Value",pch = 20)
            plot(likelihoods, main = "Red represents rejected proposals, blue represents accepted proposals",ylab= "Value",pch = 20, col = colors)
        }
        
        #function to plot mean beta values over time
        plot_beta_estimates <- function(topic){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
            }
            mean_se <- as.data.frame(matrix(0,nrow=4, ncol = 2))
            mean_se <- cbind(c("M-M", "M-F","F-M", "F-F"),mean_se)
            names(mean_se) <- c("Tie","Parameter_Estimate","SE")
            for(j in 1:length(mean_se[,1])){
                mean_se[j,2] <- mean(betas[,j])
                mean_se[j,3] <- sd(betas[,j]) 
            }
            plot <- ggplot(mean_se, aes(x=Tie, y=Parameter_Estimate))+geom_errorbar(aes(ymin=Parameter_Estimate-SE, ymax=Parameter_Estimate+SE),col = c("black","red", "green", "blue"), width=.1, size = 1.3)+geom_point(fill= c("black","red", "green", "blue"),colour="black",pch=c(21,22,23,24), size = 4)+ labs(title = paste("Topic:",topic))
            return(list(plot))
        }
        # ================================================ #
        
        pdf(paste(out_directory,county_name,"_Likelihoods.pdf", sep = ""),height=8,width=12,pointsize=7)
        par(mfrow= c(1,1))
        plot_likelihoods(Itterations) 
        dev.off()
        
        # ======= Generate Plots and save as PDFs =========#      
        make_plots <- function(Width, Height, parrow,parcol){
            #generate pdf of intercepts
            print("Plotting Intercepts...")
            pdf(file=paste(out_directory,county_name,"_Intercepts.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Topics,plot_intercepts)
            dev.off()
            
            #generate pdf of intercepts 
            print("Plotting Betas...")
            pdf(file=paste(out_directory,county_name,"_Betas.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Topics,plot_betas)
            dev.off()
            
            #generate pdf of intercepts 
            print("Plotting Latent Space Positions...")
            pdf(file=paste(out_directory,county_name,"_LS_Positions.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Topics,plot_LS_Positions)
            dev.off()
            
            #generate network plots
            print("Plotting Networks...")
            pdf(file=paste(out_directory,county_name,"_Network_Plots.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            sapply(1:Topics,plot_Topic_Network)
            dev.off()
            
            
            #generate network plots one per page
            print("Plotting Networks...")
            pdf(file=paste(out_directory,county_name,"_Network_Plots_Full_Page.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(1,1))
            sapply(1:Topics,plot_Topic_Network_Full)
            dev.off()
            
            
            #generate network plots
            print("Plotting Beta Estimates...")
            pdf(file=paste(out_directory,county_name,"_Beta_Estimates.pdf",sep = ""), width = Width, height = Height)
            par(mfrow= c(parrow,parcol))
            plots <- sapply(1:Topics,plot_beta_estimates)
            #for(j in 1:length(plots)){
            #    print(plots[[j]])
            #}
            multiplot(plotlist = plots, cols = 10)
            dev.off()
        }
        
        #stolen from the R cookbook
        multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
            require(grid)
            
            # Make a list from the ... arguments and plotlist
            plots <- c(list(...), plotlist)
            
            numPlots = length(plots)
            
            # If layout is NULL, then use 'cols' to determine layout
            if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols), byrow= T)
            }
            
            if (numPlots==1) {
                print(plots[[1]])
                
            } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                    
                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                }
            }
        }
        
        #now actually generate the plots based on the number of topics     
        if(Topics <= 25){
            make_plots(12,12,5,5)
        }else if(Topics <= 50){
            make_plots(25,12,5,10)
        }else if(Topics <= 100){
            make_plots(25,25,10,10)
        }else{
            make_plots(50,25,10,20)   
        }
        
        par(mfrow= c(1,1))
        
        
        # ======= Now generate top words and topic model diagnostics ====== #
        
        print("Generating Topic Model Output")
        
        
        
        #Establish Ordering 
        #ordering <- order(Email_Assignments, decreasing = FALSE)
        #Email_Assignments <- Email_Assignments[ordering]
        #Top_Ten_Words <- Top_Ten_Words[ordering]
        
        
        
        pdf(paste(out_directory,county_name,"_Top_Words.pdf", sep = ""),height=12,width=10,pointsize=7)
        par(mfrow= c(1,1))
        bp <- barplot2(Email_Assignments, beside = TRUE, horiz = TRUE, col = "lightblue1", border= "lightblue1", ylab = "Topic:", xlab = "Number of Edges Assigned to Topic",main = paste("Top Words by Number of Words Assigned to Topic for ",county_name,sep = "")) 
        text(0,bp,Top_Ten_Words,cex=1,pos=4)# label on the bars themselves 
        dev.off()
        
        
        
        pdf(paste(out_directory,county_name,"_Log_Ratio_LUD.pdf", sep = ""),height=8,width=12,pointsize=7)
        par(mfrow= c(2,1))
        plot_ratio_lud(Itterations) 
        dev.off()
        
        print("Outputing Geweke statistics..")
        
        intercept_list <- list()
        beta_list <- list()
        LS_list <- list()
        t= 1
        for(t in 1:Topics){
            print(paste("Current Topic: ", t))
            #Intercepts
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[Itterations+i]][t]
            }
            int <- as.numeric(geweke.diag(intercepts)$z)
            intercept_list <- append(intercept_list,int)
            
            #Betas
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][t,]
            }
            bets <- rep(0,4)
            for(i in 1:4){
                bets[i] <- as.numeric(geweke.diag(betas[,i])$z)
            } 
            beta_list <- append(beta_list,bets)
            
            #Latent Space Positions
            lat <- matrix(0,nrow=Actors, ncol= Latent_Spaces)
            for(a in 1:Actors){
                LSP <- matrix(0,nrow=Itterations, ncol= Latent_Spaces)
                for(i in 1:Itterations){
                    LSP[i,] <- Metropolis_Results[[i]][,t,a]
                }
                #calcaulte statistic
                for(j in 1:Latent_Spaces){
                    lat[a,j] <- as.numeric(geweke.diag(LSP[,j])$z)
                }
            }
            LS_list <- append(LS_list,lat)
            
        }
        
        ints <- length(which(abs(unlist(intercept_list)) < 1.645))
        print(paste("Percent Intercepts Converged ($z < 1.645$):",ints/length(intercept_list), "Number:", ints))
        
        
        bets <- length(which(abs(unlist(beta_list)) < 1.645))
        print(paste("Percent Betas Converged ($z < 1.645$):",bets/length(beta_list), "Number:", bets))
        
        ls <- length(which(abs(unlist(LS_list)) < 1.645))
        print(paste("Percent LS Positions Converged ($z < 1.645$):",ls/length(LS_list), "Number:", ls))
        
        total <- ints + bets +ls
        print(paste("Percent All Parameters Converged ($z < 1.645$):",total/(length(beta_list)+ length(LS_list) + length(intercept_list)), "Number:", total))
        
        
        #generate asortativity - token dataset
        print("Generating topic- token- assortativity dataset")
        get_beta_estimates <- function(topic){
            betas <- matrix(0,ncol = 4, nrow = Itterations)
            for(i in 1:Itterations){
                betas[i,] <- Metropolis_Results[[2*Itterations+i]][topic,]
            }
            mean_se <- as.data.frame(matrix(0,nrow=4, ncol = 2))
            mean_se <- cbind(c("M-M", "M-F","F-M", "F-F"),mean_se)
            names(mean_se) <- c("Tie","Parameter_Estimate","SE")
            for(j in 1:length(mean_se[,1])){
                mean_se[j,2] <- mean(betas[,j])
                mean_se[j,3] <- sd(betas[,j]) 
            }
            return(mean_se[,2:3])
        }
        
        beta_averages <- matrix(0, ncol = 4,nrow = Topics)
        beta_se <- matrix(0, ncol = 4,nrow = Topics)
        for(i in 1:Topics){
            result <- get_beta_estimates(i)
            beta_averages[i,] <- result[,1]
            beta_se[i,] <- result[,2]
            
        }
        
        
        topic_intercepts <- matrix(0,ncol = 2,nrow= Topics)
        for(t in 1:Topics){
            intercepts <- rep(0,Itterations)
            for(i in 1:Itterations){
                intercepts[i] <- Metropolis_Results[[Itterations+i]][t]
            }
            topic_intercepts[t,1] <- mean(intercepts)
            topic_intercepts[t,2] <- sd(intercepts)
        }
        
        topic_average_distances <- rep(0,Topics)
        for(t in 1:Topics){
            print(paste("Calculating topic latent space average distance:",t))
            #slice <- Topic_Present_Edge_Counts[,,t]
            coordinates <- matrix(0,ncol = Latent_Spaces,nrow = Actors)
            for(i in 1:Actors){
                tem<- rep(0 ,Itterations)
                for(k in 1:Itterations){
                    tem[k] <- Metropolis_Results[[k]][1,t,i]
                }
                coordinates[i,1] <- mean(tem)
                tem2<- rep(0 ,Itterations)
                for(k in 1:Itterations){
                    tem2[k] <- Metropolis_Results[[k]][2,t,i]
                }
                coordinates[i,2] <- mean(tem2)
            }
            #coordinates <- cbind(mean(Metropolis_Results[1:Itterations][1,t,]),mean(Metropolis_Results[1:Itterations][2,t,]))
            dist <- 0
            for(a in 1:Actors){
                for(b in 1:Actors){
                    if(a != b){
                        for(l in 1:Latent_Spaces){
                            dis <- (coordinates[a,1] - coordinates[b,2])^2
                        }
                        dist <- dist + sqrt(dis)
                    }
                }
            }
            dist <- dist/(Actors*(Actors -1))
            topic_average_distances[t] <- dist
        }
        
        
        
        #generate dataset for analysis
        transpose <- t(Word_Type_Topic_Counts)
        
        probs <- matrix(0, ncol = num_words, nrow = Topics)
        for(i in 1:Topics){
            probs[i,] <- Word_Type_Topic_Counts[,i]/Topic_Token_Totals[i] 
        }
        topic_data <- cbind(topic_intercepts,topic_average_distances,beta_averages,beta_se,transpose,probs)
        topic_data <- as.data.frame(topic_data)
        
        prob_vocab <- as.vector(vocab[,1])
        append_prob <- function(str){
            return(paste("PR-",str,sep = ""))
        }
        temp <- as.vector(sapply(prob_vocab,append_prob))
        vec <- c("Intercept","Intercept-SE","Average_LS_Distance","MM", "MF","FM", "FF","MM-SE", "MF-SE","FM-SE", "FF-SE",as.vector(vocab[,1]) ,temp)
        names(topic_data) <- vec
        save(topic_data,file=paste(out_directory ,county_name,"Topic_Assortativity_Data.Rdata", sep = ""))
    }
    
    
    
    
    
    
    
    
    
    
    
}#end of function definition



