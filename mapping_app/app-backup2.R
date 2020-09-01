###############################################################
######## Shiny App to display Site C Bird Distribution ########
###############################################################
{
library(shiny)
library(shinythemes)
library(tidyverse)
library(leaflet)
library(htmltools)
library(rgdal)
library(sp)
library(GISTools)
library(htmlwidgets)
library(shinyWidgets)
library(mapview)
}
  
#### Prepare Data and Base Maps ####

## Read all the data in

setwd("C:/Users/Jeffo/OneDrive/Documents/R Working/Site C Bird Mapping App")
#setwd("C:/Users/jeff.matheson/Documents/OneDrive  Personal/OneDrive/Documents/R Working/Site C Bird Mapping App")

# Songbirds

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
  # Link Songbird observations with station coordinates, rename to songdat for consistency
  songdat <- obs %>% left_join(stations, by = 'Sample Station Label') %>% 
    dplyr::select('Sample Station Label', 'Survey Name', Footprint, ValleyUpstPine, 
                  ValleyUpstDam, Code, Longitude, Latitude) %>% 
                  filter(str_detect(Code, "B-U", negate = TRUE))
  # Remove the records with no coordinates. 
  songdat <- na.omit(songdat)

  # Create the list of species from the data and exclude unknowns. 
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
  songdat <- left_join(songdat, SpeciesList, "Code")

## Woodpeckers
  woodat <- read_csv("WOODto2019.csv") 
  # Need to get lat/long coordinates to use leaflet.
  woodat <- SpatialPointsDataFrame(coords = cbind(woodat$UTM_Northing, woodat$UTM_Easting), 
                                   data = woodat, proj4string=CRS("+proj=utm +zone=10 +north +datum=NAD83"))  
  woodat <- spTransform(woodat, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  woodat <- cbind(woodat@data, woodat@coords)
  woodat <- woodat %>% rename(Longitude = coords.x1, Latitude = coords.x2)
  wooobs <- woodat %>% pivot_longer(cols = c("B-ATTW", "B-BBWO", "B-DOWO", "B-HAWO", "B-NOFL", "B-PIWO", "B-YBSA"),
                                    names_to  = 'Code', values_to = "Count") %>% 
                                    filter(Count > 0) %>% 
                                    left_join(BCspecies, by = "Code")
  wooSpList <- tibble(Code = c("B-ATTW", "B-BBWO", "B-DOWO", "B-HAWO", "B-NOFL", "B-PIWO", "B-YBSA"))
  wooSpList <- wooSpList %>% left_join(BCspecies, by = "Code") 
  
## CONI
  CONIdat <- read_csv("CONIto2019.csv") 
  CONIdat <- SpatialPointsDataFrame(coords = cbind(CONIdat$UTM_E, CONIdat$UTM_N), 
                                   data = CONIdat, proj4string=CRS("+proj=utm +zone=10 +north +datum=NAD83"))  
  CONIdat <- spTransform(CONIdat, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  CONIdat <- cbind(CONIdat@data, CONIdat@coords)
  #CONIdat <- CONIdat %>% filter()
  
## Footprint
  footprint <- readOGR("shp/fp_wgs84.shp", layer = "fp_wgs84")

# Create basemap

basemap <- leaflet() %>% 
  addTiles(group = "OSM (default)")  %>%   
  addProviderTiles(providers$Esri.WorldStreetMap, group = "Street Map") %>%    
  addProviderTiles(providers$Stamen.TonerLite, group = "Black and White") %>%
  addProviderTiles(providers$Esri.WorldImagery, group = "Imagery") %>%
  addProviderTiles(providers$Esri.WorldTopoMap, group = "Topographic") %>%
  # Layers control
  addLayersControl(
    baseGroups = c("OSM (default)", "Street Map", "Black and White", "Imagery", "Topographic"),
    overlayGroups = c("Point Count Locations","Detections"),
    options = layersControlOptions(collapsed = TRUE)) %>% 
  # Add measurement tool
  addMeasure(primaryLengthUnit = "metres",
             primaryAreaUnit = "hectares") %>% 
  # Add footprint
  addPolygons(color = "#444444", weight = 1, smoothFactor = 0.5,
              opacity = 0.8, fillOpacity = 0, data = footprint, 
              group = "Boundaries") 
basemap

#### Shiny App ####

ui <- navbarPage(
  theme = shinytheme("flatly"),
  title = "Site C Breeding Bird Follow-up Monitoring",
  id = "nav",
  collapsible = TRUE,

    tabPanel("Point Counts", 

      div(class="outer", tags$head(includeCSS("styles.css")),
                 
        leafletOutput("songmap", width = "100%", height = "100%"),

        absolutePanel(id = "controls", top = 125, left = 85, width = 250, fixed = TRUE,
                      draggable = TRUE, height = "auto",
                      h4("Breeding Bird Follow-up Monitoring Program"),
                      h4("Point Counts"),
                      p("This tool displays the location of birds detected during point count surveys
                        conducted in 2006, 2008, 2011, 2012, 2016, 2017, 2018 and 2019."),
                 
                      selectInput(inputId = "sp",
                                  label = "Select Species:",
                                  choices = SpeciesList$'English Name',
                                  selected = "Alder Flycatcher")
                      )
        )
    ),

    tabPanel("Woodpecker Surveys",
    
      div(class="outer", tags$head(includeCSS("styles.css")),
               
        leafletOutput("woomap", width = "100%", height = "100%"),
               
        absolutePanel(id = "controls", top = 125, left = 85, width = 250, fixed = TRUE,
                      draggable = TRUE, height = "auto",
                  
                      h4("Breedgin Bird Follow-up Monitoring Program"),
                      h4("Woodpeckers"),
                      p("This tool displays the location of woodpeckers detected during point count call-playback 
                      surveys conducted in 2010, 2012, 2018 and 2019."),
                  
                      selectInput(inputId = "woosp",
                                  label = "Select Species:",
                                  choices = wooSpList$'English Name',
                                  selected = "American Three-toed Woodpecker")
        )                  
      )
    ),

  tabPanel("Common Nighthawk Surveys",

    absolutePanel(id = "controls", top = 125, left = 85, width = 250, fixed = TRUE,
                   draggable = TRUE, height = "auto",
                         
    h4("Breedgin Bird Follow-up Monitoring Program"),
    h4("Common Nighthawk"),
    p("This tool displays the location of Common Nighthawk detected during surveys in 2018 and 2019."),
    p("Surveys were conducted using autonomous recording units with human listening and automated detection.")
    
    # selectInput(inputId = "sp",
    #             label = "Select Species:",
    #             choices = SpeciesList$'English Name',
    #             selected = "Alder Flycatcher"
    # )
                         
    )
  )
 
)

server <- function(input, output) {

  # Create sonbird map    
    # Reactive expresssion to subset data based on user species selection  
    spobs <- reactive({songdat[songdat$`English Name`== input$sp, ]})
  
    #spobs <- songdat[songdat$`English Name`== "Alder Flycatcher", ]  # for testing
    
    # Render basemap with stations, but observations not added yet.
    # Adding stations here forces map to propery extent and zoom level on open. 
    output$songmap <- renderLeaflet({
      basemap %>% 
        # Map station data
        addCircleMarkers(lng=~Longitude, lat=~Latitude, data=stations, 
                         popup = ~htmlEscape(stations$'Survey Name'), 
                         labelOptions(sticky=FALSE),
                         radius = 1.3, color = "blue",
                         group = "Point Count Locations",
                         stroke = TRUE)
    })

    # Map observations on basemapsong - create the reactive layer
    observe({
        leafletProxy("songmap", data = spobs()) %>% 
        clearMarkers() %>% 
        clearControls() %>% 
        # Add stations
        addCircleMarkers(lng=~Longitude, lat=~Latitude, data=stations, 
                         popup = ~htmlEscape(stations$'Survey Name'), 
                         labelOptions(sticky=FALSE),
                         radius = 1.3, color = "blue",
                         group = "Point Count Locations",
                         stroke = TRUE) %>%  
        addLegend("bottomright",
                  labels = "Point Count",
                  group = "Point Count Locations",
                  colors = "blue", 
                  opacity = 1,
                  title = "Survey Locations") %>% 
        # Add observations
        addCircleMarkers(lng=~Longitude, lat=~Latitude,
                           popup = ~htmlEscape(spobs()$'Survey Name'),
                           labelOptions(sticky=FALSE),
                           color = "#EE7621",
                           fillOpacity = 1,
                           radius = 1,
                           group = "Detections") %>% 
        addLegend("bottomright",
                    layerId = "spdetectlegend",
                    values = input$sp,
                    group = "Detections",
                    pal = colorFactor(palette = "red", domain = input$sp),
                    title = "Species",
                    opacity = 0.5)
    })

  # Create woodpecker map
    # Reactive expresssion to subset data based on user species selection
    woospobs <- reactive({wooobs[wooobs$`English Name`== input$woosp, ]})

    #woospobs <- wooobs[wooobs$`English Name`== "Hairy Woodpecker", ]
    
    # Render basemap with stations, but observations not added yet.
    # Adding stations here forces map to propery extent and zoom level on open.
    output$woomap <- renderLeaflet({
      basemap %>%
        # Map station data
        addCircleMarkers(lng=~Longitude, lat=~Latitude, data=woodat,
         #                popup = ~htmlEscape(woodat$'Station'),
                         labelOptions(sticky=FALSE),
                         radius = 1.3, color = "blue",
                         group = "Point Count Locations",
                         stroke = TRUE)
    })

    # Map observations on basemap
    observe({
      leafletProxy("woomap", data = woospobs()) %>%
        clearMarkers() %>%
        clearControls() %>%
        # Add stations
        addCircleMarkers(lng=~Longitude, lat=~Latitude, data=woodat,
        #                 popup = ~htmlEscape(stations$'Survey Name'),
                         labelOptions(sticky=FALSE),
                         radius = 1.3, color = "blue",
                         group = "Point Count Locations",
                         stroke = TRUE) %>%
        addLegend("bottomright",
                  labels = "Point Count",
                  group = "Point Count Locations",
                  colors = "blue",
                  opacity = 1,
                  title = "Survey Locations") %>%
        # Add observations
        addCircleMarkers(lng=~Longitude, lat=~Latitude,
        #                 popup = ~htmlEscape(wooobs()$'Survey Name'),
                         labelOptions(sticky=FALSE),
                         color = "#EE7621",
                         fillOpacity = 1,
                         radius = 1,
                         group = "Detections")
        addLegend("bottomright",
                  values = input$woosp,
                  group = "Detections",
                  pal = colorFactor(palette = "red", domain = input$woosp),
                  title = "Species",
                  opacity = 0.5)
    })
    
    
  # Create Common Nighthawk map
    
    
}

shinyApp(ui = ui, server = server)

