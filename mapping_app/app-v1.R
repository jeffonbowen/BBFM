###############################################################
######## Shiny App to display Site C Bird Distribution ########
###############################################################
library(shiny)
library(shinythemes)
library(tidyverse)
library(leaflet)
library(htmltools)
library(rgdal)
library(sp)
library(GISTools)
library(spdplyr)

## Read all the data in

# setwd("C:/Users/Jeffo/OneDrive/Documents/R Working/Site C Bird Mapping App")

# Stations
stations <- read_csv("Songto2019stations.csv", 
                     col_types = cols(
                       `Sample Station Photos` = col_character(),
                       `Sample Station Comments` = col_character()))

# Observations
obs <- read_csv("Songto2019observations.csv",
                col_types = cols(
                  'Visit .' = col_integer(),
                  Det_UTM_E = col_double(),
                  Det_UTM_N = col_double())) 
obs <- obs %>% dplyr::select('Sample Station Label', Species) %>% 
  rename(Code = Species)

# Link observations with station coordinates
obsloc <- obs %>% left_join(stations, by = 'Sample Station Label') %>% 
  dplyr::select('Sample Station Label', 'Survey Name', Footprint, ValleyUpstPine, 
                ValleyUpstDam, Code, Longitude, Latitude)

# Remove the records with no coordinates. 
obsloc <- na.omit(obsloc)

# Create the list of species in the data and exclude unknowns. 
SpeciesList <- as_tibble(unique(obs$Code)) %>% 
  rename(Code = value) %>% 
  filter(str_detect(Code, "B-U", negate = TRUE))
BCspecies <- read_csv("4A BC Bird List 2019-12-15.csv",
                      col_types = cols('COSEWIC' = col_character())) %>% 
  dplyr::select(ID, Code = 'Species Code', 'English Name', 'Family', 'BC List', COSEWIC, SARA)
Sort <- read_csv("4C SortSongbird.csv") %>% 
  dplyr::select(Code = 'Species Code', 'Phylo_Sort', 'IsSongbird')
SpeciesList <- SpeciesList %>% left_join(BCspecies, by = 'Code') %>% 
  left_join(Sort, 'Code')
SpeciesList <- SpeciesList %>% arrange(Code)

# Join species names
obsloc <- left_join(obsloc, SpeciesList, "Code")

# Footprint
footprint <- readOGR("shp/BB_2017_ProjectFootprint.shp", 
                     layer = "BB_2017_ProjectFootprint")
newcrs = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
footprint <- spTransform(footprint, newcrs)

# Create basemap
{
  basemap <- leaflet(options = leafletOptions(minZoom = 5, maxZoom = 30)) %>% 
    addTiles(group = "OSM (default)")  %>%   
    addProviderTiles(providers$Esri.WorldStreetMap, group = "Street Map") %>%    
    addProviderTiles(providers$Wikimedia, group = "Wikimedia") %>%
    addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") %>%
    # Layers control
    addLayersControl(
      baseGroups = c("OSM (default)", "Street Map", "Wikimedia", "Imagery"),
      overlayGroups = c("Point Count Locations","Detections"),
      options = layersControlOptions(collapsed = TRUE)) %>% 
    # Add measurement tool
    addMeasure(primaryLengthUnit = "metres",
               primaryAreaUnit = "hectares") %>% 
    # Add footprint
    addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
                opacity = 0.8, fillOpacity = 0, data = footprint, 
                group = "Boundaries") %>% 
    # Map station data
    addCircleMarkers(lng=~Longitude, lat=~Latitude, data=stations, 
                     popup = ~htmlEscape(stations$'Survey Name'), 
                     labelOptions(sticky=FALSE),
                     radius = 1.3,
                     color = "blue",
                     group = "Point Count Locations",
                     stroke = TRUE) %>%  
    addLegend("bottomright",
              labels = "Point Count",
              group = "Point Count Locations",
              colors = "blue",
              title = "Survey Locations",
              opacity = 1)
  basemap
}

## User Interface ##

ui <- fluidPage(
  
  # Set theme
  theme = shinytheme("flatly"),
  
  # App title
  titlePanel("Site C Point Count Bird Surveys 2006-2019"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    sidebarPanel(
      
      p ("This application displays the location of birds detected during point count surveys
         conducted in 2006, 2008, 2011, 2012, 2016, 2017, 2018 and 2019."),
      hr(),
      p ("Select a species below."),
      hr(),
      hr(),
      
      selectInput(inputId = "sp",
                  label = "Select Species:",
                  choices = SpeciesList$'English Name',
                  selected = "Alder Flycatcher"
      )
    ),  
    
    # Main panel for displaying outputs
    mainPanel(
      leafletOutput("mymap", width = "100%", height = "600px")
    )
  )
)


server <- function(input, output) {
  
  # Create Map    
  output$mymap <- renderLeaflet({
    
    spobs <- reactive({obsloc[obsloc$`English Name`== input$sp, ]})
    
    basemap %>% addCircleMarkers(lng=~Longitude, lat=~Latitude, data = spobs(),
                                 popup = ~htmlEscape(stations$'Survey Name'),
                                 labelOptions(sticky=FALSE),
                                 color = "#EE7621",
                                 radius = 1,
                                 group = "Detections") %>% 
      addLegend("bottomright",
                values = input$sp,
                group = "Species Detections",
                pal = colorFactor(palette = "red", domain = input$sp),
                title = "Species",
                opacity = 0.5)
  })
}

shinyApp(ui = ui, server = server)

