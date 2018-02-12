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
										 fileInput('file1', label= h4('Choose a FASTA File'))
							)
						),
						em('Choose a multiFasta file with your DNA sequences. Up to 5Mb is OK.'),style="background-color:pink;"),

						h3('Results'),
						fluidRow(
							column(3,
										 checkboxInput('withseq', 'Report sequences', T)
							),
							column(3,
										 checkboxInput('Gseq', 'Report G-sequences', F)
							)
						),
						fluidRow(
							column(3,
										 h5('Number of hits'),
										 textOutput('hits')
							),
							column(3,
										 h5('Number of Input Sequence'),
										 textOutput('Inputlength')
							)
						),

						downloadButton('downloadData', 'Download Results'),
						tableOutput('result')
	))

server = (function(input, output) {

	# Quick Score
	QuickScore <- reactive({
		resu <- c(signif(G4Hscore(toupper(gsub(' ','',input$seq0))),3),'(',length(strsplit(as.character(gsub(' ','',input$seq0)),NULL)[[1]]),')')
	})

	output$text1 <- renderText({QuickScore()})

	# Seeker part
	## Checking input text
	checkInLength <- reactive({
		checlen <- NULL
		if (nchar(input$seq)<as.numeric(input$k) & input$intype=='man') {checlen <- 'Input sequence shorter than the window size'}
		return(checlen)
	})

	checkInput <- reactive({
			chec <- grepl(paste0('[^(',paste(DNA_ALPHABET[1:15],collapse=','),')]'),gsub('[[:space:]]','',input$seq),ignore.case=T)
			if (!chec) {chectext <- 'DNA OK'} else {chectext <- 'wrong letter in your DNA'}
		return(chectext)
	})

	# importing input seq to biostring
	dataInput <- reactive({
		dataseq <- NULL
			inFile <- input$file1
			if (is.null(inFile))
			{return(NULL)}
			dataseq <- readDNAStringSet(inFile$datapath,'fasta')
		return(dataseq)
	})

	# seeking G4Hunt sequences
	dataProcessed <- reactive({
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
				}
				else
				{res <- NULL}
		}else{
			res <- NULL
		}
		return(res)
	})


	output$seqcheck <- renderText(checkInput())
	output$seqchecklen <- renderText(checkInLength())

	output$result <- renderTable(dataProcessed())
	output$Inputlength <- renderText(length(dataInput()[[1]]))
	output$hits <- renderText(length(dataProcessed()[,1]))

	output$downloadData <- downloadHandler(
		filename = function() {paste0(strsplit(basename(as.character(input$file1)),split='.',fixed=T)[[1]][1],'_hl=',input$hl,'_k=',input$k,'_G4Hseeked_',Sys.Date(),'.txt')},
		content = function(file) {
			write.table(dataProcessed(), file,sep='\t',col.names=T,row.names=F)
		})
})

# Run the application
shinyApp(ui = ui, server = server)

