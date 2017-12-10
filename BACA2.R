split.top.tables <- function(list.tables, fc = 1.5) { 
  gene.lists <- unlist(lapply(list.tables, function(d) { 
    rownames(d) <- sub("_at",'', rownames(d))
    down.genes <- rownames(d)[which(d$logFC < -log2(fc))]
    if(length(down.genes)==0) down.genes <- "1"
    up.genes <- rownames(d)[which(d$logFC >= log2(fc))]  
    if(length(up.genes)==0) up.genes <- "1"
    list(down.genes, up.genes)}), recursive = FALSE)
  # naming the gene lists 
  names(gene.lists) <- unlist(lapply(names(list.tables), function(x) paste(x,c("down","up"),sep="_")))
  print(length(gene.lists))
  print(lapply(gene.lists, length))
  gene.lists
}

DAVIDsearch <- function(gene.lists,
                        david.user,
                        david.url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/",
                        idType="AFFYMETRIX_3PRIME_IVT_ID", 
                        listType="Gene",
                        easeScore=1,
                        annotation="KEGG_PATHWAY",
                        species=NA)
{
  if(length(gene.lists) > 0 & !is.na(david.user)) {
    david.objs <- list()
    list.david.result <- list()
    for(i in 1:length(gene.lists)) {
      gene.list <- gene.lists[[i]]
      david.objs[[i]] <- DAVIDWebService$new(email=david.user, url=david.url)
      list.david.result[[i]] <- NA
      if(is.connected(david.objs[[i]])) {
        result <- addList(david.objs[[i]], gene.list, idType=idType, listName=names(gene.lists)[i], listType=listType)
        if(!is.na(species))  setCurrentSpecies(david.objs[[i]], species)
        cat("For the list ", names(gene.lists)[i], " you have:", "\n", sep="")
        cat(" - number of genes:",length(gene.list), "\n", sep="")
        cat(" - ", names(result)[1],":", result[[1]], "\n", sep="")
        cat(" - ", names(result)[2],":", length(result[[2]]), "\n", sep="")
        cat(" - species:", "\n")
        cat(getSpecieNames(david.objs[[i]]), "\n")
        setAnnotationCategories(david.objs[[i]], annotation)
        list.david.result[[i]] <- getFunctionalAnnotationChart(david.objs[[i]], threshold=easeScore, count=2)
        cat(" - number of ",annotation," terms found: ",dim(list.david.result[[i]])[1], "\n", sep="")
      }
      else
        warning(paste("Not able to connect with DAVID webservice by using:"), david.user)
    }
    return(list.david.result)
  }
  else stop("Provide a list of up-/down-regulate gene lists and a correct david user email.")
}

GOLOCsearch <- function(gene.list, 
                        annType="GO.db", 
                        annDB="org.Hs.eg.db",     
                        idType = "ENTREZID", 
                        adj.method = "BH",
                        verbose = F) 
{
  # check if the indicated packages are installed
  if(!annType %in% rownames(installed.packages()))
    stop(paste("The annotation package", annType,"is not installed",sep=" "))
  else 
    require(annType, character.only=TRUE)
  if(!annDB %in% rownames(installed.packages()))
    stop(paste("The annotation package", annDB,"is not installed",sep=" "))
  else 
    require(annDB, character.only=TRUE)
  keyType <- "GO"
  annDB <- get(annDB)
  gene.list.len <- length(gene.list)
  mig <- 1
  if (verbose)  {
    print("Valid ID types:")
    print(keytypes(annDB))
  }
  if(sum(keytypes(annDB) %in% idType) == 0) {
    stop(paste("Gene ID not valid. Available GENE IDs for the specified annotation database:", 
               paste(keytypes(org.Hs.eg.db), collapse = ","), sep="\n"))
  }
  if(mig < 0) stop("The input param mig must be >= 1.")
  if(verbose) print("Creating Res DF...")
  suppressMessages(selectDF <- select(annDB, keys=gene.list, columns=keyType, keytype=idType))
  selectDF <- selectDF[selectDF$EVIDENCE != "ND",]
  selectDF <- selectDF[,c(1,2)]
  if(length(which(is.na(selectDF[,2]))) > 0) {
    if(verbose) print("Remove NA. Unmapped to annotation...")
    selectDF <- selectDF[-which(is.na(selectDF[,2])),]
  }
  colnames(selectDF) <- c("Symbol", "ID")
  selectDF <- unique(selectDF)
  Symbol <- c()
  ID <- c()
  resDF <- plyr::ddply(selectDF, plyr::.(ID), summarise, paste(Symbol, collapse=","), length(Symbol))
  colnames(resDF) <- c("ID", "Symbols", "User_Genes")
  resDF$User_Genes_Not_In <- gene.list.len-resDF$User_Genes
  if(verbose) print("Getting All Genes for mapped IDs...")
  suppressMessages(selectAnnDF <- select(annDB, keys=resDF$ID, columns="ENTREZID", keytype=keyType))
  colnames(selectAnnDF)[1] <- "ID"
  if(verbose) print("Getting Gene count for each ID...")
  resDF$ALL_Genes <- plyr::ddply(selectAnnDF,.(ID),nrow)[,2]
  if(verbose) print("Getting Total Gene Count...")
  genome.len <- length(unique(keys(annDB, keytype="ENTREZID")))
  resDF$ALL_Genes_Not_In <- genome.len-resDF$ALL_Genes
  if(verbose) print("Calculating PValue...")
  resDF$pValue <- apply(resDF, 1, function(x){
    xx <- as.numeric(x[3:6])
    mat <- matrix(c(xx[1],xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
    testRes <- fisher.test(mat)
    testRes$p.value
  })
  if(verbose) print("Adjusting PValue...")
  resDF$adjValue <- p.adjust(resDF$pValue, method = adj.method, n = length(resDF$pValue))
  if(verbose) print("Calculating EASE...")
  resDF$ease <- apply(resDF, 1, function(x){
    xx <- as.numeric(x[3:6])
    genesIn <- xx[1]-1
    mat <- matrix(c(genesIn,xx[2],xx[3],xx[4]), nrow=2, dimnames=list(c("In", "Not_In"), c("User", "All")))
    testRes <- stats::fisher.test(mat)
    testRes$p.value
  })
  if(verbose) print("Completed. Returning...")
  resDF = resDF[-which(resDF$User_Genes <= mig),]
  rownames(resDF) = resDF[,1]
  resDF = resDF[,-1]
  more.go.info = suppressMessages(select(GO.db, keys=rownames(resDF), columns=c("TERM", "ONTOLOGY"), keytype="GOID"))
  resDF$Term = more.go.info$TERM
  resDF$Type = more.go.info$ONTOLOGY
  resDF
}

decodeStrings <- function(s, before, after, verbose=FALSE, 
                          addNames=FALSE, drop.na=TRUE, warn.if.gt.1=TRUE) 
{
  if(length(s) > 1) {
    result <- lapply(s, decodeStrings, 
                     before=before, after=after, verbose=FALSE)
    if(addNames) names(result) = s
    return(result)
  }
  starts <- (valStrings <- gregexpr(before, s)[[1]]) + attr(valStrings, "match.length")
  ends <- regexpr(after, (vStrings <- substring(s, starts,1e6)) ) - 2
  result <- substring(s, starts, starts+ends)
  if(verbose)
    cat(paste("=>",starts, starts+ends, result, sep=":", collapse="\n"))
  result <- result[result != ""]
  if(drop.na) result <- result[!is.na(result)]
  result <- result[starts>=0 & ends>=0]
  if((length(result) > 1) & (warn.if.gt.1))
    warning("More than one substring found.")
  return(result)
}
 
DAVIDquery <-function (gene.list, 
                       idType = "AFFYMETRIX_3PRIME_IVT_ID", 
                       annot = list("GOTERM_BP_FAT"), 
                       URLlengthLimit = 2048, writeHTML = TRUE,
                       verbose = FALSE) 
{
  DAVIDURLBase <- "https://david.abcc.ncifcrf.gov/"
  gene.list <- paste(gene.list, collapse = ",")
  ids <- paste(strsplit(gene.list, " ")[[1]], sep = "", collapse = "")
  firstURLOK <- FALSE
  while (firstURLOK == FALSE) {
    firstURL <- paste(DAVIDURLBase, "api.jsp?", "type=",  idType, "&ids=", ids, "&tool=", "chartReport", sep = "")
    firstURL <- paste(firstURL, "&annot=", annot, sep = "")
    if (verbose) cat("DAVIDQuery:  firstURL = ", firstURL, "\n")
    if (nchar(firstURL) < URLlengthLimit) firstURLOK <- TRUE
    else stop("The URL is too big!")
  }
  DAVIDQueryResult <- try({
    myCurlHandle <- RCurl::getCurlHandle(cookiefile = "DAVIDCookiefile.txt")
    firstStageResult <- RCurl::getURL(firstURL, curl = myCurlHandle, verbose = FALSE, ssl.verifypeer = FALSE)
    if (writeHTML)  writeChar(firstStageResult, "firstStageResult.html")
    DAVIDaction <- decodeStrings(firstStageResult, "document.apiForm.action = \"",  "\"")
    DAVIDvalues <- decodeStrings(firstStageResult, "document.apiForm.[a-z]*.value=\"", "\"", warn.if.gt.1 = FALSE)
    DAVIDfields <- decodeStrings(firstStageResult, "document.apiForm.", ".value=\"", warn.if.gt.1 = FALSE)
    secondURL <- paste(DAVIDURLBase, DAVIDaction, "?", paste(DAVIDfields, "=", DAVIDvalues, sep = "", collapse = "&"), sep = "")
    if (verbose) 
      cat("DAVIDQuery:  secondURL = ", secondURL, "\n")
    if (nchar(secondURL) > URLlengthLimit) 
      stop(paste("nchar(secondURL) too long; ", nchar(secondURL),  ">", URLlengthLimit))
    secondStageResult <- RCurl::getURL(secondURL, curl = myCurlHandle,  verbose = FALSE, ssl.verifypeer = FALSE)
    hasSessionEnded <- length(grep("Your session has ended",  secondStageResult) > 0)
    if (hasSessionEnded) 
      warning("Warning: Session ended")
    if (writeHTML) 
      writeChar(secondStageResult, "secondStageResult.html")
    downloadFileName <- decodeStrings(secondStageResult, "href=\"data/download/", "\" target=")
    if (length(downloadFileName) == 0) 
      warning("Warning: downloadFileName is not found in reply html. \n")
    downloadURL <- paste(DAVIDURLBase, "data/download/", downloadFileName, sep = "")
    if (verbose) cat("downloadURL = ", downloadURL, "\n")
    #read.delim(downloadURL, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE);
    file <- RCurl::getURL(downloadURL, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
    df <- read.delim(textConnection(file))
    for(c in 1:ncol(df)) {
      if(is.factor(df[,c]))
        df[,c] <- as.character(df[,c])
    }
    df
  })
  return(DAVIDQueryResult)
}

BBplot <- function(list.david.obj, max.pval = 0.01,
                   min.ngenes = 5, max.ngenes = 500,
                   adj.method = "Benjamini",
                   title = "BBplot",
                   name.com = "***",
                   labels = c("down", "up"),
                   colors = c("#009E73", "red"),
                   print.term = "full")
{
  if(length(list.david.obj) > 0) { # & all(unlist(lapply(list.david.obj, FUN = function(x) {inherits(x, "DAVIDFunctionalAnnotationChart")})))
    filt.dat <- lapply(list.david.obj, FUN = function(l) {
      if(adj.method == "") 
        l <- l[which(l$PValue <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "Bonferroni") 
        l <- l[which(l$Bonferroni <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "Benjamini") 
        l <- l[which(l$Benjamini <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      if(adj.method == "FDR") 
        l <- l[which(l$FDR <= max.pval & l$Count >= min.ngenes & l$Count <= max.ngenes), ]
      l <- l[which(l$PValue <= max.pval & l$Count >= min.ngenes),]
      if(print.term == "name") l$Term <- unlist(lapply(l$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
      else if(print.term == "description") l$Term <- unlist(lapply(l$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
      else if(print.term != "full") stop("Invalid input value for print.term.")
      l
    })
    list.paths <- unique(unlist(lapply(filt.dat, FUN=function(x) {unlist(x$Term)})))
    if(!is.null(list.paths)) {
      row <- rep(list.paths, length(filt.dat))
      col <- rep(1:length(filt.dat), each = length(list.paths))
      dn.up <- rep(rep(c("dn","up"), each = length(list.paths)),(length(filt.dat)/2))
      circle.size <- vector(mode = "numeric", length = (length(list.paths)*length(filt.dat)))
      significance <- vector(mode = "numeric", length = (length(list.paths)*length(filt.dat)))
      for(i in 1:length(filt.dat)) {
        sp <- ((i-1)*length(list.paths))
        for(j in 1:length(list.paths)){
          crt.path <- list.paths[j]
          ids <- which(filt.dat[[i]]$Term == crt.path)
          if(length(ids) > 0) {
            circle.size[(sp+j)] <- filt.dat[[i]]$Count[ids[1]]
            significance[(sp+j)] <- 1 - filt.dat[[i]]$PValue[ids[1]]
          }
          else {
            circle.size[(sp+j)] <- 0
            significance[(sp+j)] <- 0
          }
        }
      }
      bbplot <- NULL
      if((length(list.david.obj) %% 2) == 0 & length(colors) == 2) {
        grid.info <-  factor(rep(name.com, each = (length(list.paths)*2)), levels = name.com)
        dataset <- data.frame(row=row, col=col, circle.size=circle.size, significance=significance, 
                              dn.up=as.factor(dn.up), exp=as.factor(grid.info))
        sub.dataset <- subset(dataset, circle.size>0)
        sub.dataset$exp = droplevels(sub.dataset$exp)
        sub.dataset$dn.up = droplevels(sub.dataset$dn.up)
        if(length(levels(sub.dataset$dn.up)) == 1) {
          if(levels(sub.dataset$dn.up)[1] == "dn")
            colors = colors[1]
          else
            colors = colors[2]
        }
        print(table(sub.dataset$dn.up))
        bbplot <- ggplot(sub.dataset, aes(y = factor(row),  x = factor(col))) +
          geom_point(data = sub.dataset, aes(colour = dn.up, size = circle.size)) + #alpha=significance
          scale_colour_manual(breaks = c("dn", "up"), labels = labels, values = colors, name="Status gene") +
          scale_alpha(guide="none") +
          facet_grid(facets=.~exp, scales="free_x", space="free_x")  +
          scale_size(name = "#Genes", range = c(1, 10)) + 
          labs(title = title) +
          theme_bw() +
          theme(axis.text.x = element_blank(), plot.title = element_text(size = rel(1), colour = "blue")) +
          labs(x=NULL, y = NULL)
      }
      else {
        grid.info <- factor(rep(name.com, each = length(list.paths)), levels = name.com)
        dataset <- data.frame(row=row, col=col, circle.size=circle.size, significance=significance, exp=grid.info)
        sub.dataset <- subset(dataset, circle.size>0)
        sub.dataset$exp = droplevels(sub.dataset$exp)
        ## build the BBplot
        bbplot <- ggplot(sub.dataset, aes(y = factor(row),  x = factor(col))) +
          geom_point(data = sub.dataset, aes(size = circle.size), colour = colors[1]) + # alpha=significance
          scale_alpha(guide="none") +
          facet_grid(facets=.~exp, scales="free_x", space="free_x")  +
          scale_size(name = "#Genes", range = c(1, 10)) + 
          labs(title = title) +
          theme_bw() +
          theme(axis.text.x = element_blank(), plot.title = element_text(size = rel(1), colour = "blue")) +
          labs(x=NULL, y = NULL)
      }
      return(bbplot)
    }
    else stop("The list of selected annoations is empty.")
  }
  else stop("Provide a list of DAVIDFunctionalAnnotationChart objects.")
}

BBplotGOLOC <- function (list.enrich.res, max.pval = 0.01, 
                         min.ngenes = 5, max.ngenes = 500, 
                         adj.method = F, title = "BBplot", 
                         name.com = "***", go_category = "BP",
                         labels = c("down", "up"), 
                         colors = c("#009E73", "red"), 
                         print.term = "full") 
{
  if (length(list.enrich.res) > 0) {
    filt.dat <- lapply(list.enrich.res, FUN = function(l) {
      if (adj.method) 
        l <- l[which(l$Type == go_category & 
                       l$adjValue <= max.pval & 
                       l$User_Genes >= min.ngenes & 
                       l$User_Genes <= max.ngenes), ]
      else
        l <- l[which(l$Type == go_category & 
                       l$pValue <= max.pval & 
                       l$User_Genes >=  min.ngenes & 
                       l$User_Genes <= max.ngenes), ]
      if (print.term == "name") 
        l$Term <- rownames(l)
      else if (print.term == "id") 
        l$Term <- paste(rownames(l), l$Type, sep = ".")
      else if (print.term == "full") 
        l$Term <- paste(l$Term, l$Type, sep = ".")
      l
    })
    list.paths <- unique(unlist(lapply(filt.dat, FUN = function(x) { x$Term })))
    if (!is.null(list.paths)) {
      row <- rep(list.paths, length(filt.dat))
      col <- rep(1:length(filt.dat), each = length(list.paths))
      dn.up <- rep(rep(c("dn", "up"), each = length(list.paths)), (length(filt.dat)/2))
      circle.size <- vector(mode = "numeric", length = (length(list.paths) * length(filt.dat)))
      significance <- vector(mode = "numeric", length = (length(list.paths) *  length(filt.dat)))
      for (i in 1:length(filt.dat)) {
        sp <- ((i - 1) * length(list.paths))
        for (j in 1:length(list.paths)) {
          ids <- which(filt.dat[[i]]$Term == list.paths[j])
          if (length(ids) > 0) {
            circle.size[(sp + j)] <- filt.dat[[i]]$User_Genes[ids[1]]
            if (adj.method) 
              significance[(sp + j)] <- 1 - filt.dat[[i]]$adjValue[ids[1]]
            else
              significance[(sp + j)] <- 1 - filt.dat[[i]]$pValue[ids[1]]
          }
          else {
            circle.size[(sp + j)] <- 0
            significance[(sp + j)] <- 0
          }
        }
      }
      bbplot <- NULL
      if ((length(list.enrich.res)%%2) == 0 & length(colors) ==  2) {
        grid.info <- factor(rep(name.com, each = (length(list.paths) *  2)), levels = name.com)
        dataset <- data.frame(row = row, col = col, circle.size = circle.size, 
                              significance = significance, dn.up = dn.up, 
                              exp = grid.info)
        bbplot <- ggplot(dataset, aes(y = factor(row), x = factor(col))) + 
          geom_point(data = subset(dataset, circle.size > 0), aes(colour = dn.up, size = circle.size)) + #alpha = significance
          scale_colour_manual(breaks = c("dn",  "up"), labels = labels, values = colors, name = "Status gene") + 
          scale_alpha(guide = "none") + facet_grid(facets = . ~  exp, scales = "free_x", space = "free_x") + 
          scale_size(range = c(1, 10), name = "#genes") + labs(title = title) + 
          theme_bw() + theme(axis.text.x = element_blank(),  plot.title = element_text(size = rel(1), colour = "blue")) + 
          labs(x = NULL, y = NULL)
      }
      else {
        grid.info <- factor(rep(name.com, each = length(list.paths)),  levels = name.com)
        color.info <- rep(colors, each = length(list.paths))
        print(colors)
        dataset <- data.frame(row = row, col = col, 
                              colinfo = factor(color.info, levels = colors),
                              circle.size = circle.size, 
                              significance = significance, exp = grid.info)
        bbplot <- ggplot(dataset, aes(y = factor(row),  x = factor(col))) + 
          geom_point(data = subset(dataset, circle.size > 0), aes(size = circle.size, colour = colinfo)) + # alpha = significance
          scale_alpha(guide = "none") + facet_grid(facets = . ~ exp, scales = "free_x", space = "free_x") + 
          scale_size(range = c(1, 10), name = "#genes") + 
          scale_colour_manual(values = colors, guide = "none") + 
          labs(title = title) + theme_bw() + 
          theme(axis.text.x = element_blank(), plot.title = element_text(size = rel(1), colour = "black")) + 
          labs(x = NULL, y = NULL) 
      }
      return(bbplot)
    }
    else stop("The list of the selected annotations is empty.")
  }
  else stop("Provide a list of DAVIDFunctionalAnnotationChart objects.")
}

Jplot <- function(david.obj.1,
                  david.obj.2,
                  max.pval = 0.01,
                  min.ngenes = 5,
                  title = "Jplot",
                  print.term = "full")
{
    if(inherits(david.obj.1, "DAVIDFunctionalAnnotationChart") & inherits(david.obj.2, "DAVIDFunctionalAnnotationChart")) {
        david.obj.1 <- david.obj.1[which(david.obj.1$PValue <= max.pval & david.obj.1$Count >= min.ngenes),]
        david.obj.2 <- david.obj.2[which(david.obj.2$PValue <= max.pval & david.obj.2$Count >= min.ngenes),]
        if(print.term == "name") {
          david.obj.1$Term <- unlist(lapply(david.obj.1$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
          david.obj.2$Term <- unlist(lapply(david.obj.2$Term, FUN = function(x) unlist(strsplit(x,":"))[1]))
        }
        else if(print.term == "description") {
            david.obj.1$Term <- unlist(lapply(david.obj.1$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
            david.obj.2$Term <- unlist(lapply(david.obj.2$Term, FUN = function(x) unlist(strsplit(x,":"))[2]))
        }
        else if(print.term != "full") 
          stop("Invalid input value for cod.term. Please, indicate a number between 0 and 2.")
        size.v <- length(david.obj.1$Term)*length(david.obj.2$Term)
        row <- rep(david.obj.1$Term, each=length(david.obj.2$Term))
        col <- rep(david.obj.2$Term, length(david.obj.1$Term))
        ij <- rep(0, length(david.obj.1$Term)*length(david.obj.2$Term))
        cnt <- 1
        for(i in seq_along(david.obj.1$Term)){
            for(j in seq_along(david.obj.2$Term)) {
                genes.i <- unlist(strsplit(david.obj.1$Genes[[i]], ", "))
                genes.j <- unlist(strsplit(david.obj.2$Genes[[j]], ", "))
                ij[cnt] <- length(intersect(genes.i,genes.j))/length(union(genes.i,genes.j))
                cnt <- cnt + 1
            }
        }
        data <- data.frame(row=factor(row), col=factor(col), ij=ij)
        data$ij<-cut(data$ij,
        breaks=c(0,0.01,0.25,0.5,0.75,0.99,1),
        include.lowest=TRUE,
        label=c("0%","10-25%","25-50%","50-75%","75-99%","1"))
        jplot <- ggplot(data, aes(x=row, y=col)) +
        geom_tile(aes(fill=ij),colour="black") +
        scale_fill_brewer(palette = "YlOrRd",name="Similarity score") +
        theme(axis.text.x=element_text(angle=-90, hjust=0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text = element_text(size=20))
        return(jplot)
    }
    else stop("Provide two DAVIDFunctionalAnnotationChart objects.")
}
