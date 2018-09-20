library(shiny)
require(GenomicRanges)
library(Biostrings)

source('./seekG4hunt.r')

ui <- fluidPage(

	headerPanel("G4Hunter MultiFasta",windowTitle='Apps for G4Hunter'),
	p('by L. Lacroix, laurent.lacroix@inserm.fr'),
	helpText(a("Click Here to open the README",href="README.html",target="_blank")),

	wellPanel(h2('Quick G4Hunter calculator'),
						textInput("seq0",label= h4("Sequence"),value="GGGTTAGGGTTAGGGTTAGGG",width='100%'),
						fluidRow(
							column(2,
										 h4("Score (length)")),
							column(4,
										 h4(textOutput("text1"),style="color:red"))),
						style="background-color:lightgreen;"),

	wellPanel(style="background-color:lightblue;",h2('G4Hunter MultiFasta Seeker'),
						fluidRow(
							column(2,
										 textInput("hl",label= h4("Threshold"),value=1.5,width='80px')
							),
							column(3,
										 textInput("k",label= h4("Window size"),value=20,width='120px')
							)
						),

						wellPanel(fluidRow(
							column(4,
										 fileInput('file1', label= h4('Choose a (multi)FASTA File'))
							),

							column(3,
										 checkboxInput('withseq', 'Report sequences', T)
							),
							column(3,
										 checkboxInput('Gseq', 'Report G-sequences', F)
							)
						),
						em('Choose a multiFasta file with your DNA sequences. Up to 5Mb is OK. '),
						br(),
						em('Beware that computation can become slow if the number of entries in the multifasta file is big (>1000).'),
						h5(textOutput('Inputlength')),
						strong(textOutput('seqchecklen'),style="color: red"),
						actionButton('go','Please click here to start the computation'),
						style="background-color:pink;"),

						h3('Results'),
						h5(textOutput('hits')),

						downloadButton('downloadData', 'Download Results'),
						tableOutput('result')
	))

server = (function(input, output) {

	showModal(modalDialog(
		title = "G4Hunter notes",
		"The top part of this page allows you to compute the G4Hunter score of a single sequence. The main App (G4Hunter MultiFasta Seeker) identifies DNA regions in multiple sequences for which the G4Hunter score is above the threshold in windows of the selected size. Please cite Bedrat A, Lacroix L & Mergny JL (2016) Re-evaluation of G-quadruplex propensity with G4Hunter. Nucleic Acids Res 44(4):1746-1759, when reporting results obtained with this App.",
		easyClose = TRUE
	))

	# Quick Score
	QuickScore <- reactive({
		resu <- c(signif(G4Hscore(toupper(gsub(' ','',input$seq0))),3),'(',length(strsplit(as.character(gsub(' ','',input$seq0)),NULL)[[1]]),')')
	})

	output$text1 <- renderText({QuickScore()})

	# Seeker part
	# importing input seq to biostring
	dataInput <- reactive({
		dataseq <- NULL
			inFile <- input$file1
			if (is.null(inFile))
			{return(NULL)}
			dataseq <- readDNAStringSet(inFile$datapath,'fasta')
		return(dataseq)
	})
	## Checking number of entries
	checkInLength <- reactive({
		checlen <- NULL
		if (length(dataInput())>200) {checlen <- 'Large number of entries in your Fasta file'}
		if (length(dataInput())>1000) {checlen <- 'Very large number of entries in your Fasta file'}
		return(checlen)
	})
	output$seqchecklen <- renderText(checkInLength())

	# seeking G4Hunt sequences
	dataProcessed <- eventReactive(input$go,{
		if (!is.null(dataInput()))
		{
			#suppressWarnings()
			senam <-names(dataInput())
			hunted <- suppressWarnings(do.call(c,lapply(1:length(dataInput()),function(i) modG4huntref(k=as.numeric(input$k),hl=as.numeric(input$hl),chr=dataInput()[[i]],seqname=senam[i],with.seq=input$withseq,Gseq.only=input$Gseq))))
				if (length(hunted)!=0)
				{
					res <- as.data.frame(hunted)
					colnames(res)[8] <- 'threshold'
					colnames(res)[9] <- 'window'
					nam1 <- paste(res[,1],res[,2],sep='_')
					res <- cbind(nam1,res)
					colnames(res)[1] <- 'hitnames'
				}
				else
				{res <- NULL}
		}else{
			res <- NULL
		}
		return(res)
	})

	output$result <- renderTable(dataProcessed(),rownames=F)
	output$Inputlength <- renderText({paste0('Number of Fatsa entries: ',length(dataInput()))})
	output$hits <- renderText({paste0('Number of Hits: ',length(dataProcessed()[,1]))})

	output$downloadData <- downloadHandler(
		filename = function() {paste0(strsplit(basename(as.character(input$file1)),split='.',fixed=T)[[1]][1],'_thr=',input$hl,'_k=',input$k,'_G4Hseeked_',Sys.Date(),'.txt')},
		content = function(file) {
			write.table(dataProcessed(), file,sep='\t',col.names=T,row.names=F)
		})
})

# Run the application
shinyApp(ui = ui, server = server)

