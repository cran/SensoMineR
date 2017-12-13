format_holos <- function(path.data) {

  options(warn = -1)

  # number of subjetcs
  nbsubj <- length(list.files(path.data))

  # number of stimuli
  exemple <- read.table(paste(path.data, list.files(path.data)[1], "001_data.txt", sep = "/"), sep = ",", header = FALSE)
  exemple[, 5] <- as.factor(exemple[, 5])
  nbstim <- nlevels(exemple[, 5])

  # dataset with the names
  name.subjects <- as.data.frame(matrix(NA, 1, nbsubj))
  rownames(name.subjects) <- c("name")

  # dataset with the digit-tracking data
  raw.datadigit <- as.data.frame(matrix(NA, 1, 6))
  colnames(raw.datadigit) <- c("subject", "step", "stimulus", "time", "coordX", "coordY")
  datadigit <- list()
  steps.by.subject <- list()

  # dataset with the final configurations (coordinates)
  datafinal_coord <- as.data.frame(matrix(NA, nbstim, 2 * nbsubj))
  rownames(datafinal_coord) <- levels(exemple[, 5])
  colnames(datafinal_coord) <- paste(paste(rep("S", 2 * nbsubj), rep(1 : nbsubj, each = 2), sep = ""), rep(c("coordX", "coordY"), nbsubj), sep = "_")

  # dataset with the final configurations (verbalisation)
  datafinal_verb <- as.data.frame(matrix(NA, nbstim, nbsubj))
  rownames(datafinal_verb) <- levels(exemple[, 5])
  colnames(datafinal_verb) <- paste(rep("S", nbsubj), 1 : nbsubj, sep = "")

  # subject by subject
  for (i in 1 : nbsubj) {
    # name of the subject
    name_subj <- mixedsort(sort(list.files(path.data)))[i]
    name.subjects[1, i] <- name_subj
    colnames(name.subjects)[i] <- paste("S", i, sep = "")
    file_subj <- list.files(paste(path.data, name_subj, sep = "/"))

    # digit-tracking data
    digit_subj <- read.table(paste(path.data, name_subj, "001_data.txt", sep = "/"), sep = ",", header = FALSE)
    digit_subj <- digit_subj[, -c(1, 4)]
    colnames(digit_subj) <- c("coordX", "coordY", "stimulus", "time")
    digit_subj$stimulus <- as.factor(digit_subj$stimulus)
    step0 <- as.data.frame(matrix(NA, nlevels(digit_subj$stimulus), ncol(digit_subj)))
    for (j in 1 : nlevels(digit_subj$stimulus)) {
      stim <- levels(digit_subj$stimulus)[j]
      steps.stim <- digit_subj[which(digit_subj$stimulus == stim), ]
      step0[j, 1 : 2] <- steps.stim[1, 1 : 2]
      step0[j, 3] <- as.character(steps.stim[1, 3])
      step0[j, 4] <- as.character(steps.stim[1, 4])
    }
    datastep <- cbind.data.frame(rep(0, nrow(step0)), step0)
    colnames(datastep) <- c("step", colnames(digit_subj))
    change <- vector()
    change[1] <- 0
    stimulus_int <- as.integer(digit_subj$stimulus)
    for(j in 2 : length(stimulus_int)) {
      change[j] <- stimulus_int[j] - stimulus_int[j-1]
    }
    change[length(stimulus_int)] <- 1
    stimchange <- which(change != 0)
    step.inf <- datastep
    for(j in 1 : length(stimchange)) {
      step <- as.data.frame(matrix(NA, 1, ncol(datastep)))
      colnames(step) <- colnames(datastep)
      position <- stimchange[j]
      step$step <- j
      if(position == nrow(digit_subj)) {
        stim.step <- as.character(digit_subj[position, 3])
        step[1, 2 : 3] <- digit_subj[position, 1 : 2]
        step[1, 4] <- stim.step
        step$time <- as.character(digit_subj[position, 4])
      } else {
        stim.step <- as.character(digit_subj[position - 1, 3])
        step[1, 2 : 3] <- digit_subj[position - 1, 1 : 2]
        step[1,4] <- stim.step
        step$time <- as.character(digit_subj[position - 1, 4])
      }
      step.inf <- rbind.data.frame(step.inf, step)
    }
    datadigit_subj <- cbind.data.frame(rep(i, nrow(step.inf)), step.inf)
    colnames(datadigit_subj)[1] <- "subject"
    datadigit_subj <- datadigit_subj[, c(1, 2, 5, 6, 3, 4)]
    raw.datadigit <- rbind.data.frame(raw.datadigit, datadigit_subj)
    datadigit_subj$step <- as.factor(datadigit_subj$step)
    datadigit_subj$stimulus <- as.factor(datadigit_subj$stimulus)
    ini_step <- datadigit_subj[which(datadigit_subj$step == "0"), ]
    n_step <- datadigit_subj[which(datadigit_subj$step != "0"), ]
    steps.by.s <- as.data.frame(matrix(0, nrow(ini_step), (nrow(n_step)*2+2), list(as.character(ini_step[, "stimulus"])), byrow = TRUE))
    steps.by.s[, 1:2] <- ini_step[, c("coordX", "coordY")]
    for (j in 1 : nrow(n_step)) {
      steps.by.s[, ( 2 * j + 1) : (2 * j + 2)] <- steps.by.s[, (2 * j - 1) : (2 * j)]
      temp <- n_step[j, c("coordX", "coordY")]
      steps.by.s[as.character(n_step[j, "stimulus"]), (2 * j + 1) : (2 * j + 2)] <- temp
    }
    colnames(steps.by.s) <- paste(paste(rep("S", ncol(steps.by.s)), rep(i,ncol(steps.by.s)), sep = ""), paste(rep("step", ncol(steps.by.s) / 2 - 1), rep(0 : (ncol(steps.by.s) / 2 - 1), each = 2), rep(c("_X", "_Y"), ncol(steps.by.s) / 2 - 1), sep = ""), sep = ".")
    steps.by.subject[[i]] <- steps.by.s
    names(steps.by.subject)[i] <- paste("S", i, sep = "")

    # final configuration (coordinates)
    datadigit_subj$stimulus <- as.factor(datadigit_subj$stimulus)
    for (j in 1 : nbstim) {
      steps.stim <- datadigit_subj[which(datadigit_subj$stimulus == rownames(datafinal_coord)[j]), ]
      datafinal_coord[j, (2 * i - 1) : (2 * i)] <- steps.stim[nrow(steps.stim), 5 : 6]
    }

    # final configuration (verbalisation)
    verb_subj <- read.table(paste(path.data, name_subj, "001_comment.txt", sep = "/"), sep = ",", header = FALSE)
    verb_subj <- as.matrix(verb_subj)
    verb_subj <- verb_subj[- (1 : 6), ]
    verb_subj <- strsplit(verb_subj, "::")
    verb_subj <- unlist(verb_subj)
    list.stim <- as.character(verb_subj)
    pos.stim <- which(list.stim%in%rownames(datafinal_verb))
    change.group <- vector()
    for (j in 1 : length(pos.stim)) {
      change.group[j] <- pos.stim[j] - pos.stim[j - 1]
    }
    change.group[1] <- 0
    change.group <- pos.stim[which(change.group != 1)]
    for (j in 1 : length(change.group)) {
      list.stim.group <- vector()
      list.verb.group <- ""
      if (j != length(change.group)) {
        pos.group <- change.group[j] : (change.group[ j + 1] - 1)
      } else {
        pos.group <- change.group[j] : length(list.stim)
      }
      for (l in 1 : length(pos.group)) {
        if (list.stim[pos.group[l]]%in%rownames(datafinal_verb)) {
          list.stim.group[length(list.stim.group) + 1] <- list.stim[pos.group[l]]
        } else {
          list.verb.group <- paste(list.verb.group, verb_subj[pos.group[l]], sep = "")
        }
      }
      datafinal_verb[which(rownames(datafinal_verb)%in%list.stim.group),i] <- rep(list.verb.group, length(list.stim.group))
    }
  }
  colnames(datafinal_verb) <- paste(colnames(datafinal_verb), rep("_verb", nbsubj), sep = "")
  raw.datadigit <- raw.datadigit[-1, ]
  datadigit[[1]] <- raw.datadigit
  datadigit[[2]] <- steps.by.subject
  names(datadigit) <- c("raw.datadigit", "steps.by.subject")

  options(warn = 0)

  results <- list(name.subjects, datadigit, datafinal_coord, datafinal_verb)
  names(results) <- c("IDsubjects", "datadigit", "datafinal_coord", "datafinal_verb")

  return(results)

}
