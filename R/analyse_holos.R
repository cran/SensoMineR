analyse_holos <- function(data, method, axes = c(1, 2), graph = TRUE, export.res = FALSE) {

  results <- list()

  # functions

  plotTablet <- function(data) {

    nbsubj <- ncol(data) / 2
    for (i in 1 : nbsubj) {
       data_nappe <- data[, (i * 2 - 1) : (i * 2)]
       colnames(data_nappe) <- c("X", "Y")
	   X <- data_nappe$X
	   Y <- data_nappe$Y
       nappe <- ggplot(data_nappe, aes(x = X, y = Y, label = rownames(data_nappe))) +
           geom_text(hjust = 0, vjust = 0, colour = "blue", size = 6) +
           geom_point(colour = "blue") +
           ggtitle(paste("Final configuration - S", i, sep = "")) +
           xlim(-30, 1050) +
           ylim(-30, 710) +
           theme(
              title = element_text(size = 18, face = "bold"),
              plot.title = element_text(margin = margin(b = 20, unit = "pt")),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.background = element_blank())
        name_file <- paste("S", i, ".pdf", sep = "")
        pdf(paste(getwd(), name.folder.res, "Final configurations/Tablets", name_file, sep = "/"), width = 11.092511, height = 8.255506)
        print(nappe)
        grid.roundrect(x = 0.5, y = 0.475, width = 0.93, height = 0.87, gp = gpar(fill = "#00000000", col = "black", lwd = 70))
        grid.circle(x = 0.5, y = 0.04, r = 0.025, gp = gpar(fill = "#00000000", col = "white", lwd = 4))
        dev.off()
        if (method == "N") {
            message(paste("Exportation: ", round(i / nbsubj * 100 / 2, 1), "% done", sep = ""))
        } else {
            message(paste("Exportation: ", round(i / nbsubj * 100, 1), "% done", sep = ""))
        }
    }
  }
  summarize.digit <- function(data)  {

    nbsubj <- length(unique(data$subject))
    data$time <- as.POSIXct(data$time, format = "%H:%M:%S")
    sum1 <- data.frame(subject = numeric(), number.step = numeric(), time = numeric())
    for (i in 1 : nbsubj) {
       sum1[i, "subject"] <- i
       sum1[i, "number.step"] <- nrow(data[which(data$subject == i & data$step != "0"), ])
       j_step <- data[which(data$subject == i), ]
       dt <- difftime(as.POSIXct(j_step[nrow(j_step), "time"]),as.POSIXct(j_step[1, "time"]), units = "secs")
       sum1[i, "time"] <- format(.POSIXct(dt, tz = "GMT"), "%H:%M:%S")
    }
    sum1$subject <- as.factor(sum1$subject)
    sum1$number.step <- as.factor(sum1$number.step)
    sum2 <- list()
    for (i in 1 : nbsubj) {
       dj <- data[which(data$subject == i), ]
       freq <- as.data.frame(table(dj$stimulus))
       rownames(freq) <- freq[, 1]
       datj <- data.frame(subject = numeric(), stimulus = numeric(), freq = numeric())
       datj[1 : length(unique(dj$stimulus)), "subject"]  <- rep(i, length(unique(dj$stimulus)))
       datj[1:length(unique(dj$stimulus)), "stimulus"] <- rownames(freq)
       datj[1:length(unique(dj$stimulus)), "freq"] <- freq[2] - 1
       sum2[[i]] <- datj
    }
    results <- list(nbstep.time = sum1, freq = sum2)
    names(results$freq) <- paste(rep("S", each = nbsubj), 1 : nbsubj, sep = "")
    return(results)
 }
  plot.FA <- function(res, type, axes, cex) {

    dev.new()
    if (class(res)[[1]] == "MCA") {
        coord.group <- res$var$eta2[, axes]
        coord.ind <- res$ind$coord[, axes]
        coord.var <- res$var$coord[, axes]
		}
    if (class(res)[[1]] == "HMFA") {
        coord.group <- res$group$coord[[2]][, axes]
        coord.ind <- res$ind$coord[, axes]
        coord.var <- res$quali.var$coord[, axes]
		}
		if (class(res)[[1]] == "MFA") {
        coord.group <- res$group$coord[, axes]
        coord.ind <- res$ind$coord[, axes]
		}

    if(type == "group"){
        coord <- coord.group
        xlim <- c(0,1)
        ylim <- c(0,1)
        pch <- 17
        color <- "blue"
        title <- "Subjects representation"
    }
    if(type == "ind"){
        coord <- coord.ind
        xlim <- c((min(coord[,1]) * 1.2), (max(coord[,1]) * 1.2))
        ylim <- c((min(coord[,2]) * 1.2), (max(coord[,2]) * 1.2))
        pch <- 16
				color <- "black"
        title <- "Stimuli representation"
    }
    if(type == "var"){
        coord <- coord.var
        xlim <- c((min(coord[,1]) * 1.2),(max(coord[,1]) * 1.2))
        ylim <- c((min(coord[,2]) * 1.2),(max(coord[,2]) * 1.2))
        pch <- 16
				color <- "red"
        title <- "Stimuli representation + description"
    }
		plot(0, 0, main = title, xlab = paste("Dim ", axes[1],"(", signif(res$eig[axes[1], 2], 4), "%)", sep = ""), ylab = paste("Dim ", axes[2], " (", signif(res$eig[axes[2], 2], 4), "%)", sep = ""),xlim = xlim, ylim = ylim, col = "white", asp = 1, cex = cex)
		text(coord[, 1], y = coord[, 2], labels = rownames(coord), pos = 3, col = color, cex = cex)
		points(coord[, 1], y = coord[, 2], pch=pch, col = color, cex = cex)
		abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
		if(type == "var"){
				text(coord.ind[, 1], y = coord.ind[, 2], labels = rownames(coord.ind), pos = 3, col = "black", cex = cex)
				points(coord.ind[, 1], y = coord.ind[, 2], pch=pch, col = "black", cex = cex)
		}
		plot.save = recordPlot()
		dev.off()
		return(plot.save)
  }
  analyse.S <- function(data, axes, graph) {

    for (i in 1 : ncol(data)) {
       data[,i] <- as.factor(data[, i])
    }
    nbsubj <- ncol(data)
    res.mca <- MCA(data, axes = axes, graph=FALSE)
    plot.subj <- plot.FA(res.mca, type = "group", axes = axes, cex = 1)
    plot.stim <- plot.FA(res.mca, type = "ind", axes = axes, cex = 1)
		plot.desc <- plot.FA(res.mca, type = "var", axes = axes, cex = 0.8)
		if (export.res == TRUE) {
			name_file1 <- "Subjects.pdf"
      pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file1, sep="/"))
      print(plot.subj)
      dev.off()
      name_file2 <- "Stimuli.pdf"
      pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file2, sep="/"))
      print(plot.stim)
     	dev.off()
      name_file3 <- "Stimuli_Description.pdf"
			pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file3, sep="/"))
     	print(plot.desc)
     	dev.off()
    }
    if (graph) {
			dev.new()
      dev.new()
			print(plot.subj)
			dev.new()
			print(plot.stim)
			dev.new()
			print(plot.desc)
    }
    return(res.mca)
  }
  analyse.SN <- function(data, axes, graph) {

     nbsubj <- ncol(data) / 3
     hierar <- list(rep(c(2, 1), nbsubj), rep(2, nbsubj))
     name.group <- list(rep(c("n", "s"), nbsubj), paste(rep("S",nbsubj), 1 : nbsubj, sep = ""))
     type.group <- rep(c("c","n"), nbsubj)
     res.hmfa <- HMFA(data, H = hierar, name.group = name.group, type = type.group, axes = axes, graph = FALSE)
     plot.subj <- plot.FA(res.hmfa, type = "group", axes = axes, cex = 1)
     plot.stim <- plot.FA(res.hmfa, type = "ind", axes = axes, cex = 1)
     plot.desc <- plot.FA(res.hmfa, type = "var", axes = axes, cex = 0.8)
     if (export.res == TRUE) {
			  name_file1 <- "Subjects.pdf"
				pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file1, sep="/"))
				print(plot.subj)
				dev.off()
				name_file2 <- "Stimuli.pdf"
        pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file2, sep="/"))
	    	print(plot.stim)
     		dev.off()
     		name_file3 <- "Stimuli_Description.pdf"
        pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file3, sep="/"))
     		print(plot.desc)
			  dev.off()
     }
     if (graph) {
				dev.new()
        dev.new()
				print(plot.subj)
				dev.new()
				print(plot.stim)
				dev.new()
				print(plot.desc)
      }
      return(res.hmfa)
  }
  analyse.N <- function(data, axes, graph) {

    nbsubj <- length(data$datadigit$steps.by.subject)
    res.mfa <- MFA(data$datafinal_coord, group=rep(2,nbsubj),  type=rep("c", nbsubj), name.group=paste(rep("S",nbsubj), 1:nbsubj, sep=""), axes=axes, graph=FALSE)
    plot.subj <- plot.FA(res.mfa, type="group", axes=axes, cex=1)
		plot.stim <- plot.FA(res.mfa, type="ind", axes=axes, cex=1)
    if (export.res == TRUE) {
			name_file1 <- "Subjects.pdf"
			pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file1, sep="/"))
			print(plot.subj)
			dev.off()
			name_file2 <- "Stimuli.pdf"
      pdf(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", name_file2, sep="/"))
	    print(plot.stim)
     	dev.off()
		}
    if (graph) {
      dev.new()
			dev.new()
			print(plot.subj)
			dev.new()
			print(plot.stim)
    }
    if (export.res == TRUE) {
        for (i in 1 : nbsubj) {
           nap.panel <- data$datafinal_coord
           digit.ind <- data$datadigit$steps.by.subject[[i]]
           nap.cbind <- merge(nap.panel, digit.ind, by = "row.names")
           nom <- nap.cbind[1]
           rownames(nap.cbind) <- nom[, 1]
           nap.cbind <- nap.cbind[, -1]
           namegroup <- c(paste("S", 1:nbsubj, sep = ""), paste("Step", 0:(ncol(digit.ind)/2 - 1), sep = ""))
           group.traj = rep(2, ncol(nap.cbind)/2)
           type.traj = rep("c", ncol(nap.cbind)/2)
           ncp.traj = (nrow(nap.panel) - 1)
           num.group.sup.traj = ((ncol(nap.panel)/2) + 1):(((ncol(nap.panel)/2) + 1) + (ncol(digit.ind)/2) - 1)
           res.step <- MFA(nap.cbind, group = group.traj, type = type.traj, ncp = Inf, name.group = namegroup, axes = axes, num.group.sup = num.group.sup.traj, graph = FALSE)
           name_file <- paste("CP_","S", i, ".pdf",sep="")
           pdf(paste(getwd(), name.folder.res, "Cognitive processes", name_file, sep="/"))
           plot.MFA(res.step, choix = "group", habillage = "group", palette = palette(c(rep("transparent", (ncol(nap.panel)/2) + 1), rep("blue", ncol(digit.ind)/2))), axes = axes, title=paste("Cognitive process - S", i, sep=""))
           coord.Step <- res.step$group$coord.sup[, axes]
           lines(coord.Step, col="blue")
           dev.off()
           corr.step <- round(res.step$group$coord.sup, 4)
           corr.step <- as.data.frame(corr.step)
           newdata <- corr.step[, axes]
           newdata$id <- rownames(newdata)
           dta <- melt(newdata, id.vars="id")
           colnames(dta) <- c("step", "dim", "corr.coef")
           dta$step <- factor(dta$step, levels = rownames(corr.step))
           dtaPlot <- dta
		   corr.coef <- dtaPlot$corr.coef
           graph_TE <- ggplot(data = dtaPlot, aes(x = step, y = corr.coef, group = 1)) +
                scale_x_discrete(labels = 0:(nrow(corr.step) - 1)) +
                ggtitle(paste("Temporal evolution of the dimensions - S",i,sep="")) +
                ylim(0, 1) + geom_point(color = "blue") + geom_line(color = "blue") +
                facet_grid(dim ~ .) + theme_bw() + theme(plot.title = element_text(face="bold"))
           name_file <- paste("TE_","S", i, ".pdf",sep="")
           ggsave(name_file, path=paste(getwd(), name.folder.res, "Cognitive processes", sep="/"), width=9.44, height=7.96)
           message(paste("Exportation: ", 50+round(i/nbsubj*100/2,1), "% done", sep=""))
        }
    }
    return(res.mfa)
  }
  prepare.dataSN <- function(data_coord, data_verb) {

     nbsubj <- ncol(data_verb)
     data_coordverb <- merge(data_coord, data_verb, by="row.names")
     dataSN <- as.data.frame(matrix(NA,nrow(data_coord),nbsubj*3))
     rownames(dataSN) <- data_coordverb[,1]
     name.col <- strsplit(colnames(data_coordverb), "_")
     name.col <- unlist(name.col)
     name.col <- name.col[-which(name.col=="coordX"|name.col=="coordY"|name.col=="verb")]
     for (i in 1:nbsubj) {
        name.subj=paste("S",i,sep="")
        colnames(dataSN)[(3*i-2):(3*i)]=paste(rep(name.subj,3), c("coordX","coordY","verb"), sep="_")
        dataSN[,(3*i-2):(3*i)]=data_coordverb[,which(name.col==name.subj)]
        dataSN[,(3*i)]=as.factor( dataSN[,(3*i)])
      }
      return(dataSN)
  }

  # name of the subjetcts (correspondence email-number)
  name.subjects <- data$name.subjects
  results[[1]] <- name.subjects
  nbsubj <- ncol(name.subjects)

  # folder with the results
  if (export.res == TRUE) {
      date <- Sys.time()
      date <- gsub("[[:punct:]]", "", date)
      name.folder.res <- paste("Results_Holos", date, sep = " ")
      dir.create(paste(getwd(), name.folder.res, sep = "/"))
      dir.create(paste(getwd(), name.folder.res, "Final configurations", sep = "/"))
      dir.create(paste(getwd(), name.folder.res, "Final configurations/Tablets", sep = "/"))
      dir.create(paste(getwd(), name.folder.res, "Final configurations/Factorial analysis", sep = "/"))
  }

	# plot the final configurations on the tablet
  if (export.res == TRUE) {
      plotTablet(data$datafinal_coord)
  }

  # summarize the task for each subject
  summary.task <- summarize.digit(data$datadigit$raw.datadigit)
  results[[2]] <- summary.task

  # analysis of the final configurations
  if (method == "S") {
      res.FA <- analyse.S(data$datafinal_verb, axes = axes, graph = graph)
      results[[3]] <- res.FA
      results[[4]] <- data$datafinal_verb
  }
  if (method == "SN") {
      datafinalSN <- prepare.dataSN(data$datafinal_coord, data$datafinal_verb)
      res.FA <- analyse.SN(datafinalSN, axes = axes, graph = graph)
      results[[3]] <- res.FA
      results[[4]] <- datafinalSN
  }
  if (method == "N") {
      if (export.res == TRUE) {
          dir.create(paste(getwd(), name.folder.res, "Cognitive processes", sep = "/"))
      }
			res.FA <- analyse.N(data, axes = axes, graph = graph)
			results[[3]] <- res.FA
			results[[4]] <- list()
      results[[4]][[1]] <- data$datadigit
      results[[4]][[2]] <- data$datafinal_coord
      names(results[[4]]) <- c("digitdata","finaldata")
	}

  if (export.res == TRUE) {
    message(paste0("Results exported in: ", paste(getwd(), name.folder.res, sep = "/")))
  }

  names(results) <- c("IDsubjects", "summary.task", "res.FA", "datasets")
	return(results)
}
